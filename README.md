# Low-weight Chen iterated-integral relation model

This repository contains the **code implementation only**. The research
manuscript and its PDF are maintained separately and are intentionally not
included here.

## Program

`mini.py` constructs the shuffle, stuffle (quasi-shuffle), distribution, and
Möbius-duality relations over the fourth-root alphabet. It converts the
relations into exact SymPy linear systems, checks consistency by comparing
coefficient and augmented ranks, and solves them while leaving free variables
symbolic.

The implementation supports weights 1, 2, and 3. At weight 3, the resulting
system has one free parameter; the program does not assign it an arbitrary
value.

## Requirements

- Python 3.10 or newer
- SymPy

Install the dependency:

```bash
python -m pip install -r requirements.txt
```

## Usage

Compute through weight 2:

```bash
python mini.py 2
```

Compute through weight 3 and print each full, labelled augmented matrix:

```bash
python mini.py 3 --matrix
```

The matrix option prints the row labels, column labels, and exact matrix
`[A | b]`. Weight 3 produces a large matrix, so omit the option when only
ranks and solved values are required.

The module can also be imported. Each call to `build(n)` is independent and
automatically computes all lower weights needed by the requested system:

```python
import mini

weight_two = mini.build(2)
weight_three = mini.build(3)

print(weight_two["rank"])
print(weight_three["rank"], weight_three["augmented_rank"])
```

The returned dictionary contains `words`, `known`, `cols`, `rows`, `A`,
`b`, `rank`, and `augmented_rank`. Solved values are available in
`mini.VAL` after calling `build(n)`.

## Tests

Run the regression suite:

```bash
python -m unittest -v
```
