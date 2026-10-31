"""Host configuration, ~/.config/turbogap-mcp/hosts.toml (or TURBOGAP_MCP_HOSTS).

    [alt]
    ssh = "alt"                                   # absent: this machine
    setup = "source /l/zarrout1/toolchain/env.sh"
    bin = "/l/zarrout1/tg_builds/tg_pdf_cpu"
    runs_root = "/l/zarrout1/mcp_runs"
    test_data = "/l/zarrout1/tg_builds/turbogap_tests"
    read_roots = ["/l/zarrout1/tg_builds"]        # links may point here too
    max_ranks = 12
"""

import os
import tomllib
from dataclasses import dataclass, field
from pathlib import Path

DEFAULT_PATH = Path.home() / ".config" / "turbogap-mcp" / "hosts.toml"


@dataclass
class HostConfig:
    name: str
    bin: str
    runs_root: str
    ssh: str = ""
    ssh_command: list = field(default_factory=lambda: ["ssh"])
    setup: str = ""
    mpirun: str = "mpirun"
    test_data: str = ""
    read_roots: list = field(default_factory=list)
    max_ranks: int = 1
    timeout: float = 60.0


def config_path():
    return Path(os.environ.get("TURBOGAP_MCP_HOSTS", DEFAULT_PATH))


def load_hosts(path=None):
    path = Path(path) if path else config_path()
    if not path.is_file():
        return {}
    raw = tomllib.loads(path.read_text())
    hosts = {}
    for name, table in raw.items():
        unknown = set(table) - set(HostConfig.__dataclass_fields__)
        if unknown:
            raise ValueError(f"{path}: [{name}] has unknown key(s) {', '.join(sorted(unknown))}")
        missing = {"bin", "runs_root"} - set(table)
        if missing:
            raise ValueError(f"{path}: [{name}] needs {', '.join(sorted(missing))}")
        hosts[name] = HostConfig(name=name, **table)
    return hosts
