"""The keyword table against the Fortran it is parsed from, read a second way."""

import re

ALL_BRANCHES_RE = re.compile(r"keyword\s*==\s*['\"](\w+)['\"]")
SUBROUTINE_RE = re.compile(r"subroutine (read_options_\w+)\b(.*?)end subroutine", re.S)


def fortran_input_keywords(tree):
    """Every spelling a read_options_* branch tests, with continuation lines joined
    and comment lines dropped: an independent reading of what TurboGAP accepts."""
    text = (tree / "src" / "read_files.f90").read_text()
    code = "\n".join(line for line in text.split("\n") if not line.lstrip().startswith("!"))
    code = re.sub(r"&[ \t]*\n[ \t]*&?", "", code)
    names = set()
    for _, body in SUBROUTINE_RE.findall(code):
        names.update(ALL_BRANCHES_RE.findall(body))
    return names


def test_every_fortran_keyword_is_known(table, tree):
    missing = sorted(n for n in fortran_input_keywords(tree) if not table.is_known(n))
    assert missing == []


def test_no_keyword_the_fortran_lacks(table, tree):
    fortran = fortran_input_keywords(tree) | set(table.prepass)
    extra = sorted(n for n in table.by_name if n not in fortran)
    assert extra == []


def test_prepass_keywords(table):
    assert table.prepass == ["pot_file", "n_species"]
    assert table.lookup("pot_file")[0]["type"] == "string"


def test_counts_match_the_reference_header(table, tree):
    md = (tree / "docs" / "keywords.md").read_text()
    n_input = int(re.search(r"\*\*The input file\*\* &mdash; (\d+) keywords", md).group(1))
    n_gap = int(re.search(r"\*\*The potential \(\.gap\) file\*\* &mdash; (\d+) keywords", md).group(1))
    assert len(table.input_keywords) == n_input
    assert len(table.gap_keywords) == n_gap


def test_lookup_matches_the_markdown_row(table, tree):
    md = (tree / "docs" / "keywords.md").read_text()
    row = next(line for line in md.splitlines() if line.startswith("| `pair_distribution_kde_sigma` |"))
    cells = [c.strip() for c in row.strip("|").split("|")]
    found = table.lookup("pair_distribution_kde_sigma")[0]
    assert found["type"] == cells[1]
    assert f"`{found['default']}`" == cells[2]
    assert found["units"] == cells[3]


def test_alias_resolves_to_its_keyword(table):
    assert table.lookup("exp_energy_scales_beg")[0]["name"] == "exp_energy_scales"
    assert table.lookup("input_file")[0]["name"] == "atoms_file"


def test_every_list_has_a_count(table):
    lists = [k.name for k in table.input_keywords if "list" in k.type or "(" in k.type]
    # box_scaling_factor takes 1, 3 or 9 values and reads them itself.
    uncounted = [n for n in lists if n not in table.counts and n != "box_scaling_factor"]
    assert uncounted == []
    assert table.counts["mc_planes"] == "4*params%mc_n_planes"
    assert table.counts["exp_labels"] == "params%n_exp"


def test_gap_lookup_has_one_entry_per_block(table):
    groups = {e["group"] for e in table.lookup("rcut", file="gap")}
    assert groups == {"soap_turbo", "distance_2b", "angle_3b"}


def test_search_ranks_name_matches_and_filters_by_mode(table):
    found = table.search("kde")
    assert found["keywords"][0]["name"] == "pair_distribution_kde_sigma"
    md_only = [k.name for k in table.input_keywords if k.modes == ["md"]]
    predict = {k["name"] for k in table.search("", mode="predict", limit=1000)["keywords"]}
    assert md_only and not predict & set(md_only)
