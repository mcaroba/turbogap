"""Where TurboGAP runs. Every operation is a bash script fed to `bash -s`; a host
only decides how that reaches the machine (here, or over one reused ssh
connection), so jobs never see ssh. Every call has a timeout: an unreachable host
is an error, not a hang.
"""

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
        return ["bash", "-s"]

    def sh(self, script, timeout=None):
        """Run a bash script on the host and return its stdout. The script is one
        group, read whole before it runs, so a command reading stdin (mpirun) cannot
        swallow the lines after it."""
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

    def reachable(self):
        try:
            self.sh("true", timeout=15)
            return True, ""
        except HostError as e:
            return False, str(e)


class SshHost(Host):
    def argv(self):
        c = self.config
        return [*c.ssh_command,
                "-o", "BatchMode=yes",
                "-o", "ConnectTimeout=8",
                "-o", "ControlMaster=auto",
                "-o", "ControlPath=~/.ssh/cm-turbogap-mcp-%r@%h:%p",
                "-o", "ControlPersist=10m",
                c.ssh, "bash -s"]


def make_host(config: HostConfig):
    return SshHost(config) if config.ssh else Host(config)
