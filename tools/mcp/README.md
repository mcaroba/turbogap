# turbogap-mcp

An MCP server for using TurboGAP: keyword lookup, input validation and runs on
configured hosts. The live calculator comes later.

| tool | does |
| --- | --- |
| `keywords` | describe one keyword (`name=`) or search them (`query=`, `mode=`, `file="input"\|"gap"`) |
| `validate_input` | check an input file for a mode, as `src/read_files.f90` will read it |
| `run_start` | stage and launch a run on a host, from an input or a regression `case`; returns a job id |
| `run_status` | running / done / failed / cancelled / lost, exit code, log tail |
| `run_results` | one trajectory frame and thermo.log, summarised; files stay on the host |
| `run_cancel` | stop a job's whole process group |

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

The second line adds real runs, which must match the regression suite's own run of
each case byte for byte; without the two variables they skip.
