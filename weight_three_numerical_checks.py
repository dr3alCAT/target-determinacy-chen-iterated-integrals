#!/usr/bin/env python3
"""Numerical spot checks for the weight-three target integrals.

This script is an independent computational aid. It checks the seven
imaginary-part relations and compares numerical quadrature with the closed
forms. These numerical checks are not part of the mathematical proof.

Dependency: mpmath (install with: python -m pip install mpmath).
"""

from mpmath import (
    mp, mpf, mpc, quad, log, pi, atan, catalan, polylog, nstr
)


def I3(a, b, c):
    """Evaluate the word integral I(a,b,c) by nested quadrature."""
    inner = lambda t: quad(lambda s: log(1 - s / a) / (s - b), [0, t])
    return quad(lambda t: inner(t) / (t - c), [0, mpf(1) / 2, 1])


def main():
    mp.dps = 30
    L, P, G = log(mpf(2)), pi, catalan
    ii, m1, mi = mpc(0, 1), mpc(-1), mpc(0, -1)
    z = mpc(0)

    W = {
        1: (ii, ii, z),
        2: (ii, z, mi),
        3: (ii, ii, m1),
        4: (ii, ii, mi),
        5: (ii, mi, mi),
        6: (ii, mi, ii),
        7: (m1, ii, ii),
        8: (ii, m1, ii),
        9: (m1, ii, m1),
        10: (ii, m1, m1),
    }
    v = {k: I3(*w).imag for k, w in W.items()}

    # Left and right sides of the seven imaginary-part relations.
    relations = [
        (-v[2], -G * L / 2 + P**3 / 192),
        (v[1] + v[2] - v[3] - v[4] - v[5], P * L**2 / 8),
        (v[6] - 2 * v[5], -G * L / 2 + 7 * P * L**2 / 32 + P**3 / 128),
        (v[7] + v[3] - v[9] - v[10], mpf(0)),
        (v[8] + 2 * v[7], 3 * P * L**2 / 32 + P**3 / 384),
        (v[7] + v[8] + v[3], P * L**2 / 8),
        (v[5] - v[6] - v[4], -P * L**2 / 32 - P**3 / 128),
    ]
    print("Seven imaginary-part relations:")
    for k, (lhs, rhs) in enumerate(relations, 1):
        print(
            f"R{k}",
            nstr(lhs, 14),
            nstr(rhs, 14),
            nstr(abs(lhs - rhs), 3),
        )

    C_value = quad(
        lambda t: t * log(1 + t) ** 2 / (1 + t**2), [0, 1]
    )
    B_value = quad(
        lambda t: atan(t)
        * (log(1 + t**2) / t - log(1 + t) / (1 + t)),
        [0, 1],
    )
    q1 = quad(
        lambda t: atan(t) * log(1 + t**2) / t, [0, 1]
    ) / 2  # Im I(i,i,0)
    q2 = quad(
        lambda t: log(1 + t) ** 2 / (1 + t**2), [0, 1]
    ) / 2  # Im I(-1,-1,i)
    s3 = polylog(3, (1 + ii) / 2).imag

    print("\nTarget-integral checks:")
    print(
        "Ccal",
        nstr(C_value, 20),
        nstr(abs(C_value - (L**3 / 24 + P**2 * L / 96)), 3),
    )
    print(
        "Bcal",
        nstr(B_value, 20),
        nstr(abs(B_value - (P**3 / 128 - 3 * P * L**2 / 32)), 3),
    )
    print(
        "Bcal word form",
        nstr(abs(2 * v[1] - v[9] - v[10] - B_value), 3),
    )
    print(
        "(iv)",
        nstr(abs(2 * q1 + q2 - (P**3 / 128 + P * L**2 / 32)), 3),
    )
    print(
        "q1",
        nstr(abs(q1 - (s3 + G * L / 2 - P * L**2 / 32 - 3 * P**3 / 128)), 3),
    )
    print(
        "q2",
        nstr(abs(q2 - (-2 * s3 - G * L + 3 * P * L**2 / 32 + 7 * P**3 / 128)), 3),
    )


if __name__ == "__main__":
    main()
