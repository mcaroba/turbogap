"""Where TurboGAP runs. Every operation is a bash script fed to `bash -s`; a host
only decides how that reaches the machine (here, or over one reused ssh
connection), so jobs never see ssh. Every call has a timeout: an unreachable host
is an error, not a hang. A host with a [slurm] table submits instead of launching.
"""

import datetime
import os
import re
import shlex
import subprocess

from .config import HostConfig


class HostError(RuntimeError):
    pass


class Host:
    def __init__(self, config: HostConfig):
        self.config = config
        self.name = config.name

    def argv(self):
        c = self.config
        if not c.ssh:
            return ["bash", "-s"]
        return [*c.ssh_command,
                "-o", "BatchMode=yes",
                "-o", "ConnectTimeout=8",
                "-o", "ControlMaster=auto",
                "-o", "ControlPath=~/.ssh/cm-turbogap-mcp-%r@%h:%p",
                "-o", "ControlPersist=10m",
                c.ssh, "bash -s"]

    def check_certificate(self):
        """An expired ssh certificate fails as a bare 'Permission denied', so it is
        caught here, before connecting, and named."""
        c = self.config
        if not c.certificate:
            return
        path = os.path.expanduser(c.certificate)
        try:
            shown = subprocess.run(["ssh-keygen", "-L", "-f", path], capture_output=True, text=True,
                                   timeout=10).stdout
        except (OSError, subprocess.TimeoutExpired):
            return
        found = re.search(r"Valid: from \S+ to (\S+)", shown)
        if not found:
            return
        until = datetime.datetime.fromisoformat(found.group(1))
        if until <= datetime.datetime.now():
            hint = f": renew it with `{c.renew}`" if c.renew else ""
            raise HostError(f"{self.name}: the ssh certificate {c.certificate} expired at {until}{hint}")

    def sh(self, script, timeout=None):
        """Run a bash script on the host and return its stdout. The script is one
        group, read whole before it runs, so a command reading stdin (mpirun) cannot
        swallow the lines after it."""
        self.check_certificate()
        wrapped = "{\nset -e\n" + script + "\n} < /dev/null\n"
        try:
            done = subprocess.run(self.argv(), input=wrapped, capture_output=True,
                                  text=True, timeout=timeout or self.config.timeout)
        except subprocess.TimeoutExpired as e:
            raise HostError(f"{self.name}: no answer within {e.timeout:.0f} s") from None
        if done.returncode != 0:
            raise HostError(f"{self.name}: exit {done.returncode}: {done.stderr.strip()[-500:]}")
        return done.stdout

    def put(self, text, path):
        """Write text to a file on the host, through a quoted heredoc whose
        delimiter the text does not contain, so nothing in it is interpreted."""
        delimiter = "TURBOGAP_MCP_EOF"
        while delimiter in text:
            delimiter += "_"
        body = text if text.endswith("\n") else text + "\n"
        self.sh(f"cat > {shlex.quote(path)} <<'{delimiter}'\n{body}{delimiter}\n")

    def exists(self, paths):
        """{path: bool} for every path, in one round trip."""
        if not paths:
            return {}
        script = "".join(f"[ -e {shlex.quote(p)} ] && echo 1 || echo 0\n" for p in paths)
        return dict(zip(paths, (f == "1" for f in self.sh(script).split())))

    def launch(self, job_dir, mode, ranks):
        """Start TurboGAP detached in a session of its own. Its pid, which is also
        its process group, goes to job_dir/pid and its exit code to job_dir/rc."""
        c = self.config
        run = f"{shlex.quote(c.bin)} {shlex.quote(mode)}"
        if ranks > 1:
            run = f"{shlex.quote(c.mpirun)} -np {int(ranks)} {run}"
        inner = "\n".join(filter(None, [
            "echo $$ > pid",
            c.setup,
            f"export OMP_NUM_THREADS={int(c.omp_threads)}",
            f"{run} > out.log 2>&1",
            "echo $? > rc",
        ]))
        self.sh(f"cd {shlex.quote(job_dir)}\n"
                f"setsid nohup bash -c {shlex.quote(inner)} > /dev/null 2>&1 < /dev/null &\n"
                "for i in $(seq 50); do [ -s pid ] && break; sleep 0.1; done\n")

    def poll(self, job_dir):
        """(state, exit code). State is running, done, cancelled, cancelling, or
        lost: gone without writing an exit code. rc is looked at again after the
        liveness test, which the job can pass, write rc and exit in between."""
        words = self.sh(f"""cd {shlex.quote(job_dir)}
c=0; [ -f cancelled ] && c=1
if [ -f rc ]; then echo done $c $(cat rc)
elif [ -s pid ] && kill -0 $(cat pid) 2>/dev/null; then echo running $c
elif [ -f rc ]; then echo done $c $(cat rc)
else echo lost $c
fi
""").split()
        state, cancelled = words[0], words[1] == "1"
        code = int(words[2]) if len(words) > 2 else None
        if cancelled:
            state = "cancelling" if state == "running" else "cancelled"
        return state, code

    def tail(self, path, lines):
        q = shlex.quote(path)
        return self.sh(f"if [ -f {q} ]; then tail -n {int(lines)} {q}; fi\n")

    def kill(self, job_dir):
        """TERM the job's process group, then KILL whatever is left after 5 s."""
        self.sh(f"""cd {shlex.quote(job_dir)}
touch cancelled
[ -s pid ] || exit 0
g=$(cat pid)
kill -TERM -- -$g 2>/dev/null || exit 0
for i in $(seq 50); do kill -0 -- -$g 2>/dev/null || exit 0; sleep 0.1; done
kill -KILL -- -$g 2>/dev/null || true
""")

    def describe(self, job_dir):
        """Scheduler details for run_status; none here."""
        return {}

    def reachable(self):
        try:
            self.sh("true", timeout=15)
            return True, ""
        except HostError as e:
            return False, str(e)


# squeue states before a job has a node, and sacct states that end one without its
# own exit code reaching rc.
QUEUED = ("PENDING", "CONFIGURING", "REQUEUED", "RESV_DEL_HOLD", "REQUEUE_HOLD", "SUSPENDED")
FAILED = ("TIMEOUT", "FAILED", "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE")


class SlurmHost(Host):
    """Jobs go through sbatch. The batch script writes rc as a launched job does;
    squeue says queued or running, and sacct says how a job without rc ended."""

    def job_script(self, job_dir, mode, ranks):
        c = self.config
        s = c.slurm
        lines = ["#!/bin/bash",
                 f"#SBATCH --job-name=tgmcp-{os.path.basename(job_dir)}",
                 f"#SBATCH --account={s['account']}",
                 f"#SBATCH --partition={s['partition']}",
                 f"#SBATCH --time={s['time']}",
                 f"#SBATCH --ntasks={int(ranks)}",
                 f"#SBATCH --cpus-per-task={int(c.omp_threads)}",
                 f"#SBATCH --output={job_dir}/out.log"]
        lines += [f"#SBATCH {option}" for option in s.get("options", [])]
        lines += [f"cd {shlex.quote(job_dir)}",
                  c.setup,
                  f"export OMP_NUM_THREADS={int(c.omp_threads)}",
                  f"srun {shlex.quote(c.bin)} {shlex.quote(mode)}",
                  "echo $? > rc"]
        return "\n".join(filter(None, lines)) + "\n"

    def scheduler(self, script, timeout=None):
        """A script that calls Slurm, after setup: some sites put Slurm in a module."""
        return self.sh(f"{self.config.setup}\n{script}", timeout)

    def launch(self, job_dir, mode, ranks):
        self.put(self.job_script(job_dir, mode, ranks), f"{job_dir}/job.sh")
        self.scheduler(f"cd {shlex.quote(job_dir)}\nsbatch --parsable job.sh | cut -d';' -f1 > slurm_id\n"
                "[ -s slurm_id ]\n")

    def poll(self, job_dir):
        """(state, exit code), as Host.poll, plus queued: submitted and waiting for a node."""
        out = self.scheduler(f"""cd {shlex.quote(job_dir)}
c=0; [ -f cancelled ] && c=1
id=$(cat slurm_id)
if [ -f rc ]; then echo "done $c $(cat rc)"; exit 0; fi
s=$(squeue -h -j "$id" -o %T 2>/dev/null || true)
if [ -n "$s" ]; then echo "slurm $c $s"; exit 0; fi
if [ -f rc ]; then echo "done $c $(cat rc)"; exit 0; fi
a=$(sacct -n -X -P -j "$id" -o State,ExitCode 2>/dev/null | head -1 || true)
echo "ended $c ${{a:-UNKNOWN|}}"
""").strip().splitlines()[-1]
        kind, cancelled, rest = out.split(" ", 2)
        cancelled = cancelled == "1"
        if kind == "done":
            return ("cancelled" if cancelled else "done"), int(rest)
        if kind == "slurm":
            if cancelled:
                return "cancelling", None
            return ("queued" if rest.split()[0] in QUEUED else "running"), None
        state, _, exit_code = rest.partition("|")
        state = state.split()[0] if state.strip() else "UNKNOWN"
        code = int(exit_code.split(":")[0]) if exit_code[:1].isdigit() else None
        if cancelled or state == "CANCELLED":
            return "cancelled", code
        if state == "COMPLETED":
            return "done", code
        if state in FAILED:
            return "failed", code
        return "lost", code

    def describe(self, job_dir):
        out = self.scheduler(f"""cd {shlex.quote(job_dir)}
id=$(cat slurm_id)
echo TURBOGAP_MCP_SLURM
echo "$id"
s=$(squeue -h -j "$id" -o '%T|%r|%N|%M' 2>/dev/null || true)
if [ -n "$s" ]; then echo "$s"; else
sacct -n -X -P -j "$id" -o State,ExitCode,NodeList,Elapsed 2>/dev/null | head -1 || true; fi
""").split("TURBOGAP_MCP_SLURM\n")[-1].splitlines()
        info = {"slurm_id": out[0]}
        if len(out) > 1:
            fields = out[1].split("|")
            if len(fields) == 4 and ":" not in fields[1]:
                info.update(slurm_state=fields[0], reason=fields[1], nodes=fields[2], elapsed=fields[3])
            elif len(fields) == 4:
                info.update(slurm_state=fields[0], slurm_exit=fields[1], nodes=fields[2], elapsed=fields[3])
        return info

    def kill(self, job_dir):
        """scancel, then wait up to 10 s for the queue to let go of the job."""
        self.scheduler(f"""cd {shlex.quote(job_dir)}
touch cancelled
id=$(cat slurm_id)
scancel "$id" 2>/dev/null || true
for i in $(seq 50); do [ -z "$(squeue -h -j "$id" 2>/dev/null || true)" ] && exit 0; sleep 0.2; done
""", timeout=self.config.timeout + 15)


def make_host(config: HostConfig):
    return SlurmHost(config) if config.slurm else Host(config)
