"""Runs: staging a job directory on a host, launching it, and the records that let
a job outlive the server. One JSON file per job, under TURBOGAP_MCP_STATE or
~/.local/state/turbogap-mcp/jobs.
"""

import json
import os
import re
import secrets
import shlex
import time
from dataclasses import asdict, dataclass
from pathlib import Path, PurePosixPath

from . import examples
from .hosts import Host
from .results import frame_script, summarise_frame, summarise_thermo
from .validate import as_dicts, parse_entries, validate

# What the regression suite adds to every staged case so its outputs do not
# depend on the node's memory (tests/regression/run.sh, stage()).
CASE_SUFFIX = "\nmax_Gbytes_per_process = 1.0\n"
FILE_KEYWORDS = ("atoms_file", "pot_file", "exp_data_files", "mc_molecule_files")
NAME_RE = re.compile(r"^[\w.+-]+$")


class JobError(ValueError):
    pass


@dataclass
class Job:
    id: str
    host: str
    dir: str
    mode: str
    ranks: int
    created: float
    case: str = ""


def state_dir():
    path = Path(os.environ.get("TURBOGAP_MCP_STATE",
                               Path.home() / ".local" / "state" / "turbogap-mcp")) / "jobs"
    path.mkdir(parents=True, exist_ok=True)
    return path


def save(job):
    (state_dir() / f"{job.id}.json").write_text(json.dumps(asdict(job), indent=1))


def load(job_id):
    if not NAME_RE.match(job_id):
        raise JobError(f"no job '{job_id}'")
    path = state_dir() / f"{job_id}.json"
    if not path.is_file():
        raise JobError(f"no job '{job_id}'")
    return Job(**json.loads(path.read_text()))


def all_jobs():
    jobs = [Job(**json.loads(p.read_text())) for p in state_dir().glob("*.json")]
    return sorted(jobs, key=lambda j: j.created)


def new_id():
    return time.strftime("%Y%m%d-%H%M%S-") + secrets.token_hex(2)


def inside(path, roots):
    p = PurePosixPath(path)
    if not p.is_absolute() or ".." in p.parts:
        return False
    return any(p == PurePosixPath(r) or PurePosixPath(r) in p.parents for r in roots if r)


def referenced_files(input_text):
    """Every file the input names, relative to the job directory unless absolute."""
    entries, _ = parse_entries(input_text)
    names = []
    for entry in entries:
        if entry.keyword in FILE_KEYWORDS:
            names.extend(v for v in entry.values if v[:1] not in ("!", "#") and v != "none")
    return names


def atoms_file_name(input_text):
    entries, _ = parse_entries(input_text)
    names = [e.values[0] for e in entries if e.keyword in ("atoms_file", "input_file") and e.values]
    if not names:
        raise JobError("atoms_xyz was given but the input names no atoms_file")
    return names[-1]


def plan_links(host: Host, links, case_data):
    """[(source on host, name in job dir)], every source inside a readable root."""
    c = host.config
    roots = [c.runs_root, c.test_data, *c.read_roots]
    out = []
    for link in links:
        source = link["source"]
        target = link.get("target") or PurePosixPath(source).name
        if case_data and not source.startswith("/"):
            source = f"{c.test_data}/{case_data}/{source}"
        if not inside(source, roots):
            raise JobError(f"link source {source} is outside {host.name}'s readable roots {roots}")
        if not NAME_RE.match(target):
            raise JobError(f"link target '{target}' must be a plain file name")
        out.append((source, target))
    return out


def stage_case(root, case):
    recipe = examples.read_case(root, case)
    files = {}
    extra = root / "tests" / "regression" / "cases" / case / "files"
    if extra.is_dir():
        files = {p.name: p.read_text() for p in sorted(extra.iterdir()) if p.is_file()}
    return recipe, recipe["input"] + CASE_SUFFIX, files


def start(host: Host, table, root, mode="", input_text="", atoms_xyz="", links=(), files=None,
          ranks=0, case=""):
    files = dict(files or {})
    case_data = ""
    if case:
        recipe, case_input, case_files = stage_case(root, case)
        mode = mode or recipe["mode"]
        input_text = input_text or case_input
        links = list(links) or recipe["links"]
        ranks = ranks or recipe["ranks"]
        files = {**case_files, **files}
        case_data = recipe["data"]
        if not host.config.test_data:
            raise JobError(f"{host.name} has no test_data configured, so case data cannot be linked")
    ranks = ranks or 1
    if not mode or not input_text:
        raise JobError("give mode and input_text, or a case")
    if ranks > host.config.max_ranks:
        raise JobError(f"{ranks} ranks is over {host.name}'s max_ranks = {host.config.max_ranks}")
    issues = validate(table, input_text, mode)
    errors = [i for i in issues if i.severity == "error"]
    if errors:
        return {"started": False, "issues": as_dicts(errors)}

    if atoms_xyz:
        files.setdefault(atoms_file_name(input_text), atoms_xyz)
    for name in files:
        if not NAME_RE.match(name):
            raise JobError(f"file name '{name}' must be a plain file name")
    planned = plan_links(host, links, case_data)
    found = host.exists([source for source, _ in planned])
    missing = [s for s, ok in found.items() if not ok]
    if missing:
        raise JobError(f"missing on {host.name}: {', '.join(missing)}")

    job_id = new_id()
    job = Job(job_id, host.name, f"{host.config.runs_root}/{job_id}", mode, ranks, time.time(), case)
    script = [f"mkdir -p {shlex.quote(job.dir)}", f"cd {shlex.quote(job.dir)}"]
    script += [f"ln -sfn {shlex.quote(s)} {shlex.quote(t)}" for s, t in planned]
    host.sh("\n".join(script) + "\n")
    host.put(input_text, f"{job.dir}/input")
    for name, text in files.items():
        host.put(text, f"{job.dir}/{name}")

    needed = [n if n.startswith("/") else f"{job.dir}/{n}" for n in referenced_files(input_text)]
    absent = [p for p, ok in host.exists(needed).items() if not ok]
    if absent:
        raise JobError(f"the input names files that are not there: {', '.join(absent)} "
                       f"(job directory {job.dir} left for inspection)")
    host.launch(job.dir, mode, ranks)
    save(job)
    warnings = [i for i in issues if i.severity == "warning"]
    return {"started": True, "job_id": job.id, "host": host.name, "dir": job.dir,
            "mode": mode, "ranks": ranks, "warnings": as_dicts(warnings)}


def status(host: Host, job: Job):
    state, code = host.poll(job.dir)
    thermo = host.sh(f"if [ -f {shlex.quote(job.dir)}/thermo.log ]; then "
                     f"grep -vc '^#' {shlex.quote(job.dir)}/thermo.log || true; fi\n").strip()
    out = {"job_id": job.id, "host": job.host, "dir": job.dir, "mode": job.mode,
           "state": state, "exit_code": code, "elapsed_s": round(time.time() - job.created),
           "out_log_tail": host.tail(f"{job.dir}/out.log", 20)}
    if thermo:
        out["thermo_rows"] = int(thermo)
    if state == "done" and code != 0:
        out["state"] = "failed"
    return out


def results(host: Host, job: Job, frame=-1, arrays=False):
    trajectory = f"{job.dir}/trajectory_out.xyz"
    thermo = f"{job.dir}/thermo.log"
    out = {"job_id": job.id, "files": {"trajectory": trajectory, "thermo": thermo,
                                       "log": f"{job.dir}/out.log"}}
    present = host.exists([trajectory, thermo])
    if present[trajectory]:
        text = host.sh(frame_script(trajectory, frame), timeout=300)
        walked, _, body = text.partition("\n")
        out["frames_walked"] = int(walked or 0)
        out["frame"] = summarise_frame(body, arrays)
    if present[thermo]:
        out["thermo"] = summarise_thermo(host.sh(f"cat {shlex.quote(thermo)}\n", timeout=300))
    if not any(present.values()):
        out["note"] = "no trajectory_out.xyz or thermo.log yet"
    return out


def cancel(host: Host, job: Job):
    host.kill(job.dir)
    state, code = host.poll(job.dir)
    return {"job_id": job.id, "state": state, "exit_code": code}
