"""The MCP surface: tool, resource and prompt registration only."""

import atexit
import functools
from functools import cache
from typing import Any, Literal

from mcp.server.mcpserver import MCPServer
from mcp.server.mcpserver.exceptions import ToolError

from . import examples, jobs
from .calculator import Registry
from .config import config_path, load_hosts
from .hosts import HostError, make_host
from .keywords import KeywordTable
from .tree import root
from .validate import as_dicts, validate

INSTRUCTIONS = """\
TurboGAP is a molecular-dynamics, Monte Carlo and prediction engine for SOAP-based
Gaussian Approximation Potentials. An input file is one `keyword = value` per line;
an unknown keyword aborts the run. Before writing an input, look keywords up with
`keywords`, start from a worked example under turbogap://examples, and check the
result with `validate_input` for the mode it will run in (predict, md, mc, soap, ipi).
Runs are asynchronous: `run_start` returns a job id at once; poll it with `run_status`
and read a summary with `run_results`. turbogap://hosts lists where runs can go.
For many single-point evaluations of one system, `calc_start` keeps TurboGAP running
with the potential loaded and `calc_evaluate` costs one force call per geometry.
"""

MAD_OBSERVABLES = {
    "xps": "xps",
    "xrd": "xrd",
    "pdf": "pair_distribution",
    "saxs": "saxs",
    "structure_factor": "structure_factor",
}

server = MCPServer("turbogap", instructions=INSTRUCTIONS)
calculators = Registry()
atexit.register(calculators.stop_all)


@cache
def table():
    return KeywordTable(root())


@server.tool()
def keywords(query: str = "", name: str = "", mode: str = "",
             file: Literal["input", "gap"] = "input", limit: int = 25) -> dict[str, Any]:
    """Look up TurboGAP keywords.

    With `name`, describe that keyword exactly (type, default, units, modes, what it
    needs, how many values a list takes, full description); a .gap keyword can match
    one entry per descriptor block. Otherwise search names and descriptions for every
    word of `query`, optionally only keywords that apply to `mode`. `file` picks the
    input file or the potential (.gap) file.
    """
    if name:
        found = table().lookup(name, file)
        if not found:
            return {"found": [], "did_you_mean": table().close_names(name.lower())}
        return {"found": found}
    return table().search(query, mode, file, limit)


@server.tool()
def validate_input(input_text: str, mode: str) -> dict[str, Any]:
    """Check a TurboGAP input file, as read_files.f90 will read it, for `mode`.

    Errors are what aborts the run or silently misreads it: an unknown keyword, a
    keyword glued to its '=', a value of the wrong type, a list shorter than its
    count (Fortran keeps reading into the next line), an unquoted '/' cutting a path
    short, a missing atoms_file or pot_file. Warnings are a keyword documented for
    another mode, one whose prerequisite is unset, a deprecated or ignored one, and
    a duplicate. `issues` is empty when the input is clean.
    """
    issues = validate(table(), input_text, mode)
    errors = sum(i.severity == "error" for i in issues)
    return {"ok": errors == 0, "errors": errors, "warnings": len(issues) - errors,
            "issues": as_dicts(issues)}


def reported(tool):
    """Pass expected failures to the client as their own message; the SDK
    otherwise keeps any exception's text on the server."""
    @functools.wraps(tool)
    def wrapper(*args, **kwargs):
        try:
            return tool(*args, **kwargs)
        except (HostError, ValueError) as e:
            raise ToolError(str(e)) from e
    return wrapper


def get_host(name):
    hosts = load_hosts()
    if name not in hosts:
        raise ValueError(f"no host '{name}' in {config_path()}; configured: {', '.join(hosts) or 'none'}")
    return make_host(hosts[name])


def job_and_host(job_id):
    job = jobs.load(job_id)
    return job, get_host(job.host)


@server.tool()
@reported
def run_start(host: str, mode: str = "", input_text: str = "", case: str = "",
              atoms_xyz: str = "", links: list[dict[str, str]] | None = None,
              files: dict[str, str] | None = None, ranks: int = 0) -> dict[str, Any]:
    """Start TurboGAP on a configured host and return a job id at once.

    Either give `mode` and `input_text`, or name a regression `case` (see
    turbogap://examples) to run it exactly as the suite stages it; anything else
    given overrides the case. `atoms_xyz` is extended-XYZ text written to the input's
    atoms_file; `files` are further small text files by name; `links` are
    {source, target} symlinks to files already on the host, inside its runs_root,
    test_data or read_roots. The input is validated first and nothing starts on an
    error; every file the input names must exist before the launch. `ranks` is
    capped by the host's max_ranks.
    """
    return jobs.start(get_host(host), table(), root(), mode=mode, input_text=input_text,
                      atoms_xyz=atoms_xyz, links=links or [], files=files, ranks=ranks, case=case)


@server.tool()
@reported
def run_status(job_id: str) -> dict[str, Any]:
    """State of a job (running, done, failed, cancelled, lost), its exit code, the
    last 20 lines of its log, and how many thermo.log rows it has written."""
    job, h = job_and_host(job_id)
    return jobs.status(h, job)


@server.tool()
@reported
def run_results(job_id: str, frame: int = -1, arrays: bool = False) -> dict[str, Any]:
    """Summary of a job's output: one trajectory frame (0-based, -1 the last) with
    energy, virial, lattice, species and force statistics, and thermo.log as decimated
    series. `arrays` adds positions and forces for up to 1000 atoms. Full files stay
    on the host; their paths are returned."""
    job, h = job_and_host(job_id)
    return jobs.results(h, job, frame, arrays)


@server.tool()
@reported
def run_cancel(job_id: str) -> dict[str, Any]:
    """Stop a job: TERM its process group, then KILL what is left after 5 s."""
    job, h = job_and_host(job_id)
    return jobs.cancel(h, job)


@server.tool()
@reported
def calc_start(host: str, input_text: str, atoms_xyz: str, links: list[dict[str, str]] | None = None,
               files: dict[str, str] | None = None, ranks: int = 1,
               start_timeout_s: float = 600.0, eval_timeout_s: float = 600.0) -> dict[str, Any]:
    """Start a live calculator: `turbogap ipi` on a host, kept running with its
    potential loaded, answering one geometry at a time over an i-PI socket.

    `input_text` is an ipi-mode input (pot_file, n_species, species, masses, e0,
    atoms_file, ...) without ipi_address, which the calculator sets. `atoms_xyz` is
    written to its atoms_file and fixes the atoms: every later geometry must have the
    same species in the same order. `links` and `files` are as for run_start. Returns
    once TurboGAP has loaded the potential and taken the handshake, which includes one
    discarded evaluation of atoms_xyz. The calculator lives as long as this server.
    """
    return calculators.start(get_host(host), table(), root(), input_text, atoms_xyz,
                             links=links or [], files=files, ranks=ranks,
                             start_timeout=start_timeout_s, eval_timeout=eval_timeout_s)


@server.tool()
@reported
def calc_evaluate(calc_id: str, xyz: str, arrays: bool = True) -> dict[str, Any]:
    """Energy (eV), forces (eV/A) and virial (eV, TurboGAP's convention) for one
    extended-XYZ frame (Lattice, species, pos) on a running calculator. The cell may
    change between calls; the atoms may not. Forces come back for up to 1000 atoms,
    statistics always. Positions travel in bohr, so they are evaluated as
    (x / 0.5291772109) * 0.5291772109, which can differ from x in the last bit."""
    return calculators.get(calc_id).evaluate(xyz, arrays)


@server.tool()
@reported
def calc_stop(calc_id: str) -> dict[str, Any]:
    """Stop a calculator: send EXIT, wait for TurboGAP to finish, close the socket."""
    return calculators.stop(calc_id)


@server.resource("turbogap://hosts", mime_type="application/json")
def host_list() -> list:
    """The configured hosts, and whether each answers now."""
    out = []
    for name, config in load_hosts().items():
        ok, why = make_host(config).reachable()
        out.append({"name": name, "ssh": config.ssh or "(this machine)", "bin": config.bin,
                    "runs_root": config.runs_root, "max_ranks": config.max_ranks,
                    "omp_threads": config.omp_threads,
                    "reachable": ok, "error": why})
    return out


@server.resource("turbogap://jobs", mime_type="application/json")
def job_list() -> list:
    """Every job this server has started, oldest first."""
    return [vars(j) for j in jobs.all_jobs()]


@server.resource("turbogap://keywords.md", mime_type="text/markdown")
def keywords_md() -> str:
    """The full keyword reference, docs/keywords.md."""
    return (root() / "docs" / "keywords.md").read_text()


@server.resource("turbogap://examples", mime_type="application/json")
def example_list() -> list:
    """Every regression case, with its mode and rank count: worked example inputs."""
    return examples.list_cases(root())


@server.resource("turbogap://examples/{case}", mime_type="application/json")
def example(case: str) -> dict:
    """One regression case: its input, mode, ranks and the data files it links in."""
    return examples.read_case(root(), case)


@server.prompt()
def mad_setup(observable: str) -> str:
    """Build a Molecular Augmented Dynamics input driven by an experimental
    observable: xps, xrd, pdf, saxs or structure_factor."""
    label = MAD_OBSERVABLES.get(observable.lower(), observable.lower())
    cases = [c["name"] for c in examples.list_cases(root())
             if observable.lower() in c["name"] and c["mode"] == "md"]
    start = ", ".join(f"turbogap://examples/{c}" for c in cases) or "turbogap://examples"
    return f"""\
Write a TurboGAP md input that biases the dynamics towards an experimental {observable} spectrum.

1. Read a working example first: {start}.
2. Call `keywords` with query "{label}" and mode "md", then with query "exp", to see every
   keyword this observable uses. The MAD block is n_exp, exp_labels (here '{label}'),
   exp_data_files, exp_n_samples, exp_energy_scales, exp_energies and exp_forces; every
   exp_* list takes n_exp values.
3. Ask me for the experimental data file, the structure (atoms_file), the potential
   (pot_file) and species, masses and e0, unless I have already given them.
4. Call `validate_input` with mode "md" and fix every error before showing me the input.
"""


def main():
    server.run()


if __name__ == "__main__":
    main()
