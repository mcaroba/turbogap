# turbogap-mcp

An MCP server for using TurboGAP. Phase 1: keyword lookup and input validation;
runs and the live calculator come later.

| tool | does |
| --- | --- |
| `keywords` | describe one keyword (`name=`) or search them (`query=`, `mode=`, `file="input"\|"gap"`) |
| `validate_input` | check an input file for a mode, as `src/read_files.f90` will read it |

| resource | is |
| --- | --- |
| `turbogap://keywords.md` | `docs/keywords.md` |
| `turbogap://examples` | every regression case, with its mode |
| `turbogap://examples/{case}` | one case: input, mode, ranks and the data it links |

Prompt `mad_setup(observable)` walks through a MAD input for xps, xrd, pdf, saxs or
structure_factor.

The keyword table is `tools/keyword_docs.py`, imported, so the server and
`docs/keywords.md` cannot disagree. The tree is the one this directory sits in,
or `TURBOGAP_ROOT`.

```sh
cd tools/mcp
uv sync --extra test
uv run pytest                                   # needs gfortran for the tokenizer oracle
claude mcp add turbogap -- "$PWD/.venv/bin/turbogap-mcp"
```
