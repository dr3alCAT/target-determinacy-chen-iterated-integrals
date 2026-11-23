# Target Determinacy in Weight-Two Chen Iterated Integrals

**Finite relation systems and exact algebraic reduction**

The manuscript source and compiled PDF are prepared locally; the repository currently contains this project description. The manuscript files will be added in a follow-up push.

## Abstract

We develop an exact relation-based method for evaluating Chen iterated integrals with letters in the fourth-root alphabet `{0} ∪ μ₄`. The construction combines shuffle and stuffle identities with distribution relations and Möbius duality, producing a finite linear system with exact algebraic coefficients and right-hand sides built from lower-weight values. The central question is target determinacy: whether a particular linear combination is determined by the available relations even when the system has a nontrivial nullspace. The weight-two calculation establishes the complete reduction of the sixteen convergent words in the alphabet. Two weight-three challenges are included at the end as open extensions.

## Files

- `main.tex` — editable LaTeX source of the latest edited manuscript.
- `main.pdf` — compiled manuscript (15 pages).

## Build

Requires a LaTeX distribution with the packages used in the source, including `amsmath`, `amssymb`, `mathtools`, `geometry`, `microtype`, `fancyhdr`, `tcolorbox`, `enumitem`, `graphicx`, `xcolor`, and TikZ.

Compile with:

```bash
latexmk -pdf main.tex
```

## Status

The weight-two relation system and its reductions are the results developed in the manuscript. The two weight-three integrals at the end are posed as challenges; no evaluations are claimed there.

## Author

Toluwanimi Charles
