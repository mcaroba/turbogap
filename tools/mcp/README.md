# turbogap-mcp

An MCP server for using TurboGAP: keyword lookup, input validation, runs on
configured hosts, and a live calculator.

| tool | does |
| --- | --- |
| `keywords` | describe one keyword (`name=`) or search them (`query=`, `mode=`, `file="input"\|"gap"`) |
| `validate_input` | check an input file for a mode, as `src/read_files.f90` will read it |
| `run_start` | stage and launch a run on a host, from an input or a regression `case`; returns a job id |
| `run_status` | running / done / failed / cancelled / lost, exit code, log tail |
| `run_results` | one trajectory frame and thermo.log, summarised; files stay on the host |
| `run_cancel` | stop a job's whole process group |
| `calc_start` | `turbogap ipi` on a host, potential loaded, behind an i-PI socket; returns a calc id |
| `calc_evaluate` | energy, forces and virial for one extxyz frame of the same atoms |
| `calc_stop` | EXIT, and wait for TurboGAP to finish |

| resource | is |
| --- | --- |
| `turbogap://keywords.md` | `docs/keywords.md` |
| `turbogap://examples` | every regression case, with its mode |
| `turbogap://examples/{case}` | one case: input, mode, ranks and the data it links |
| `turbogap://hosts` | configured hosts and whether each answers |
| `turbogap://jobs` | every job started, oldest first |

Prompt `mad_setup(observable)` walks through a MAD input for xps, xrd, pdf, saxs or
structure_factor.

The keyword table is `tools/keyword_docs.py`, imported, so the server and
`docs/keywords.md` cannot disagree. The tree is the one this directory sits in,
or `TURBOGAP_ROOT`.

## Hosts

`~/.config/turbogap-mcp/hosts.toml` (or `TURBOGAP_MCP_HOSTS`), one table per host:

```toml
[alt]
ssh = "alt"                                   # leave out for this machine
setup = "source /l/zarrout1/toolchain/env.sh"
bin = "/l/zarrout1/tg_builds/tg_pdf_cpu"
runs_root = "/l/zarrout1/mcp_runs"            # one directory per job
test_data = "/l/zarrout1/tg_builds/turbogap_tests"
read_roots = ["/l/zarrout1/tg_builds"]        # links may point here too
max_ranks = 12
omp_threads = 1                               # per rank; results depend on it
```

A host with a `[name.slurm]` table runs jobs through `sbatch` (`srun` inside), and
`run_status` reads `squeue` and `sacct`: a job is `queued` until it has a node, and
one Slurm ended (TIMEOUT, OUT_OF_MEMORY, NODE_FAIL) is `failed` with Slurm's state.
`run_cancel` is `scancel`. Slurm commands run after `setup`.

```toml
[roihu]
ssh = "roihuc"
certificate = "~/.ssh/id_lumi_new_ed25519-cert.pub"
renew = "roihucer"                            # named when the certificate expires
setup = "source /usr/share/lmod/lmod/init/bash; ..."
bin = "/scratch/project_2003566/zarrout1/turbogap_cpu_pdf"
runs_root = "/scratch/project_2003566/zarrout1/mcp_runs"
max_ranks = 384

[roihu.slurm]
account = "project_2003566"
partition = "small"
time = "00:30:00"
options = ["--mem-per-cpu=2G"]                # optional further #SBATCH lines
```

With `certificate`, an expired one is reported as such before ssh is tried.
Calculators need the socket on the node TurboGAP runs on, so they are refused on
Slurm hosts.

A calculator is a job in mode ipi whose socket is held by this server, so it ends
with the server: TurboGAP aborts when the socket closes. The server side of the
protocol is `tests/ipi_socket/refserver.py`, imported. On a remote host `ssh -R`
forwards `/tmp/ipi_<name>` there to the listener here; the host needs nothing open.

A job runs detached (`setsid`), so it outlives the server; its record is in
`~/.local/state/turbogap-mcp/jobs` (or `TURBOGAP_MCP_STATE`). Over ssh one
connection is reused and every call times out, so a host that is down is an error.

## Tests

```sh
cd tools/mcp
uv sync --extra test
uv run pytest            # gfortran for the tokenizer oracle
TURBOGAP_MCP_TEST_BIN=/path/to/turbogap TURBOGAP_DATA_ROOT=/path/to/turbogap_tests uv run pytest
claude mcp add -s user turbogap -- "$PWD/.venv/bin/turbogap-mcp"
```

The second line adds real runs: jobs must match the regression suite's own run of
each case byte for byte, and the calculator must agree with `turbogap predict` on
the geometries it sent, as `tests/ipi_socket/run.sh` checks `turbogap ipi`. The
binary needs the ipi mode. Without the two variables they skip.
