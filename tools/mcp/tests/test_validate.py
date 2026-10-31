"""The validator against inputs TurboGAP is known to accept, against mutations of
them that it rejects or misreads, and its tokenizer against gfortran itself."""

import shutil
import subprocess

import pytest

from turbogap_mcp import examples
from turbogap_mcp.validate import expand_repeats, split_list_directed, validate


def errors(issues):
    return [i for i in issues if i.severity == "error"]


def messages(issues, keyword):
    return [i.message for i in issues if i.keyword == keyword]


def test_every_regression_input_is_clean(table, tree):
    """The suite runs every one of these, so none may carry an error."""
    failed = {}
    for case in examples.list_cases(tree):
        recipe = examples.read_case(tree, case["name"])
        found = errors(validate(table, recipe["input"], recipe["mode"]))
        if found:
            failed[case["name"]] = found
    assert failed == {}


@pytest.fixture
def co_predict(tree):
    return examples.read_case(tree, "co_predict")["input"]


@pytest.mark.parametrize("line, keyword, fragment", [
    ("md_stpe = 1.0", "md_stpe", "did you mean md_step"),
    ("md_step=1.0", "md_step=1.0", "one token"),
    ("print_progress = yes", "print_progress", "not a valid logical"),
    ("write_xyz = 1.5", "write_xyz", "not a valid integer"),
    ("md_step 1.0", "md_step", "expected '='"),
    ("md_step =", "md_step", "no value"),
])
def test_mutation_is_an_error(table, co_predict, line, keyword, fragment):
    issues = errors(validate(table, co_predict + line + "\n", "predict"))
    assert any(fragment in m for m in messages(issues, keyword)), issues


def test_short_list_is_an_error(table, co_predict):
    text = co_predict.replace("e0 = -.16138053 0.", "e0 = -.16138053")
    assert any("needs 2 value(s)" in m for m in messages(errors(validate(table, text, "predict")), "e0"))


def test_list_count_follows_its_keyword(table, tree):
    text = examples.read_case(tree, "mc_planes")["input"]
    assert errors(validate(table, text, "mc")) == []
    short = text.replace("mc_planes = 0.0 0.0 1.0 8.6", "mc_planes = 0.0 0.0 1.0")
    assert any("needs 4 value(s)" in m for m in messages(validate(table, short, "mc"), "mc_planes"))


def test_unquoted_slash_is_an_error(table, co_predict):
    text = co_predict.replace('pot_file = "gap_files/CO.gap"', "pot_file = gap_files/CO.gap")
    found = messages(errors(validate(table, text, "predict")), "pot_file")
    assert any("read as 'gap_files'" in m for m in found), found


def test_upper_case_prepass_keyword_is_an_error(table, co_predict):
    text = co_predict.replace("n_species = 2", "N_SPECIES = 2")
    assert messages(errors(validate(table, text, "predict")), "N_SPECIES")


def test_missing_required_is_an_error(table, co_predict):
    text = co_predict.replace('atoms_file = "atoms.xyz"\n', "")
    assert messages(errors(validate(table, text, "predict")), "atoms_file")


@pytest.mark.parametrize("line, keyword, fragment", [
    ("md_nsteps = 5", "md_nsteps", "documented for mode(s) md"),
    ("xps_sigma = 0.2", "xps_sigma", "needs do_xps"),
    ("random_seed = 7", "random_seed", "also set on line"),
])
def test_mutation_is_a_warning(table, co_predict, line, keyword, fragment):
    issues = validate(table, co_predict + line + "\n", "predict")
    assert errors(issues) == []
    assert any(fragment in m for m in messages(issues, keyword)), issues


def test_inline_comment_is_not_a_value(table, co_predict):
    text = co_predict + "thermostat = berendsen ! the usual\nmd_step = 1.0 ! fs\n"
    assert errors(validate(table, text, "md")) == []
    assert not messages(validate(table, text, "md"), "thermostat")


def test_unknown_mode(table, co_predict):
    assert errors(validate(table, co_predict, "relax"))


RECORDS = [
    "kw = 'a b' c,d",
    'kw = "it""s" x',
    "kw = 1 2 / 3",
    "kw = gap_files/CO.gap",
    "kw = .true. F t",
    "kw = 3*1.5 2",
    "kw=3 4",
]

PROBE = """\
program probe
   character(len=256) :: record
   character(len=64) :: words(8)
   character(len=1024) :: line
   integer :: i, n, ios
   open (10, file='records.txt')
   open (11, file='counts.txt')
   do
      read (11, *, iostat=ios) n
      if (ios /= 0) exit
      read (10, '(A)') record
      words = '@'
      read (record, *, iostat=ios) words(1:n)
      line = ''
      do i = 1, n
         line = trim(line)//trim(words(i))//'|'
      end do
      write (*, '(A)') trim(line)
   end do
end program probe
"""


@pytest.mark.skipif(shutil.which("gfortran") is None, reason="gfortran not found: tokenizer oracle not run")
def test_tokenizer_matches_gfortran(tmp_path):
    """gfortran reads each record list-directed as character tokens, one more than
    the Python tokenizer finds; the tokens must agree and the extra one must stay
    unread ('@')."""
    (tmp_path / "probe.f90").write_text(PROBE)
    subprocess.run(["gfortran", "probe.f90", "-o", "probe"], cwd=tmp_path, check=True)
    (tmp_path / "records.txt").write_text("\n".join(RECORDS) + "\n")
    expected = []
    counts = []
    for record in RECORDS:
        tokens = expand_repeats(split_list_directed(record)[0])
        counts.append(str(len(tokens) + 1))
        expected.append("|".join(tokens + ["@"]) + "|")
    (tmp_path / "counts.txt").write_text("\n".join(counts) + "\n")
    got = subprocess.run(["./probe"], cwd=tmp_path, check=True, capture_output=True, text=True)
    assert got.stdout.splitlines() == expected
