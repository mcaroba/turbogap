"""The MCP surface: tool, resource and prompt registration only."""

from functools import cache
from typing import Any, Literal

from mcp.server.mcpserver import MCPServer

from . import examples
from .keywords import KeywordTable
from .tree import root
from .validate import as_dicts, validate

INSTRUCTIONS = """\
TurboGAP is a molecular-dynamics, Monte Carlo and prediction engine for SOAP-based
Gaussian Approximation Potentials. An input file is one `keyword = value` per line;
an unknown keyword aborts the run. Before writing an input, look keywords up with
`keywords`, start from a worked example under turbogap://examples, and check the
result with `validate_input` for the mode it will run in (predict, md, mc, soap, ipi).
"""

MAD_OBSERVABLES = {
    "xps": "xps",
    "xrd": "xrd",
    "pdf": "pair_distribution",
    "saxs": "saxs",
    "structure_factor": "structure_factor",
}

server = MCPServer("turbogap", instructions=INSTRUCTIONS)


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
