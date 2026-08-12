# Documenting a keyword

`docs/keywords.md`, `docs/keywords.html` and `turbogap --help` are all
generated from `src/read_files.f90` by `tools/keyword_docs.py`. Nothing in
them is written by hand, so there is one place to edit and no table to keep
in step.

Two kinds of keyword are covered, and the same `!>` block documents both:

- **input-file** keywords, parsed by the `read_options_*` family and grouped
  by which of those subroutines claims them;
- **potential-file** keywords, parsed by `read_gap_hypers` and grouped by the
  `.gap` block they appear in — `soap_turbo`, `distance_2b`, `angle_3b`,
  `core_pot`. The same name means different things in different blocks
  (`rcut`, `delta`, `sigma`, `Z1`), so these are keyed by block, and the
  defaults come from that block's derived type rather than from
  `input_parameters`.

## Adding a keyword

Put a `!>` block immediately above the branch that handles it:

```fortran
      !> @kw mc_move_max
      !> Largest displacement of a single-atom "move". The trial displacement
      !> is drawn uniformly up to this.
      !> @units A
      !> @modes mc
      !> @needs mc_types
      else if (keyword == 'mc_move_max') then
```

then `make docs`. `make check-docs` fails if a keyword has no block or if one
of the three outputs is stale; it also runs from `.pre-commit-config.yaml`, so
a commit that touches `src/read_files.f90` or `src/types.f90` regenerates them.

## Tags

| Tag | Meaning |
| --- | --- |
| `@kw` | Required, first line. Every spelling the branch tests, comma separated, in the order the source tests them. The generator compares this against the branch and **aborts** if they disagree, which is what stops a block outliving the keyword it describes. |
| `@units` | Physical units, free text. Omit for dimensionless and logical keywords. |
| `@modes` | Which of `predict`, `md`, `mc`, `soap` the keyword does anything in. Omit when it applies to all of them. Drives the filtering in `turbogap --help <mode>`. Rejected on a `.gap` keyword — the potential file is read the same way in every mode, and `turbogap --help gap` selects it instead. |
| `@needs` | Other keywords that must be set for this one to have an effect. |
| `@see` | Related keywords. |
| `@status` | `ignored` for a keyword still parsed for compatibility whose value nothing consults, `deprecated` for one superseded but still working. |
| `@type`, `@default` | Overrides, for the handful of keywords whose field cannot be found automatically. Do not use these otherwise — see below. |

## What not to write

The generator already reads these out of the source, and repeating them in a
comment is how a table starts disagreeing with the code:

- the **group**, from the enclosing `read_options_*` subroutine;
- the **fields** the branch assigns, from `params%<field>`;
- the **type** and the **default**, from the field's declaration in
  `src/types.f90`;
- the **side effects** — every other field the branch writes, reported as
  "sets X", naming the keyword that owns X where there is one. This is how
  `do_xrd` documents that it also switches on `do_pair_distribution` and
  `do_structure_factor` without anyone having to remember to say so.

## The `.gap` file

`read_gap_hypers` is two nested loops rather than a flat `else if` chain: a
`gap_beg <descriptor>` dispatcher, and inside each arm a
`do while (keyword /= "gap_end")` over that block's keywords. The generator
follows exactly that structure, so a keyword belongs to the block whose loop
encloses it. `soap_turbo` has two such loops — a pre-scan that finds
`n_species`, because none of the per-species arrays can be allocated without
it, and then the real loop — and a keyword found in both is attributed to the
first.

`gap_beg` and `gap_end` are delimiters rather than settings and are described
in the section preamble instead of being listed as keywords.

One thing the tables cannot tell you, because it is a property of the parser:
an unrecognised keyword in the **input** file aborts the run, while an
unrecognised keyword in a **`.gap`** block is skipped without comment. A
misspelling in a potential file therefore shows up as a wrong answer, not an
error.

## Why `read_files.f90` and not `types.f90`

The obvious alternative is to document each field where it is declared. It
does not work, because the fields and the keywords are not the same set:

- of the 251 fields in `input_parameters`, 239 branches parse keywords, and
  the overlap is partial in both directions — a third of the fields are
  internal state no keyword reaches, and a keyword can write several fields;
- in the `.gap` file the mismatch is worse still: `rcut`, `delta` and `sigma`
  each name a different field in a different type depending on the block they
  sit in, so there is no single declaration a comment could hang off;
- a keyword can be spelled more than one way (`atoms_file` / `input_file`,
  `verbosity` / `verb`), and the aliases exist only in `read_files.f90`;
- a keyword need not be named after the field it writes — `nd_n_samples`
  writes `xrd_n_samples`, which a comment on `xrd_n_samples` could not say;
- the side effects, the validation, the deprecations and the mode gating are
  all in `read_files.f90`.

So the prose goes where the keyword is, and the type and the default are
harvested from `types.f90`. Each fact stays in the one file that owns it, and
the parts that can drift are the parts nobody writes twice.

## Keywords marked `ignored`

Several keywords are parsed, warned about and discarded: their field was
removed from `input_parameters` because nothing read it. Deleting the keyword
as well is not an option — `read_input_file` treats an unrecognised keyword as
fatal, so every deck that sets one would abort. `ignored_keyword` consumes the
line and says so, and the reference lists them with an **accepted and ignored**
tag rather than pretending they work.
