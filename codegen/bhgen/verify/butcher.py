"""Exact verification of Runge-Kutta order conditions.

Everything here is done in :class:`fractions.Fraction`, so a "pass" is an exact
arithmetic identity, not a numerical near-miss.  This is the cheapest possible
formal check in the whole project and it is the one that catches the most
damaging class of bug: a single mistyped tableau digit silently costs several
orders of accuracy while the method still looks like it is working.

The conditions are generated from rooted trees rather than copied from a table,
so the same code verifies any explicit tableau to any order.  For a method of
order p one needs, for every rooted tree t with at most p nodes,

    sum_i b_i Psi_i(t) = 1 / gamma(t)

where Psi_i(bullet) = 1, Psi_i(t) = prod over subtrees u of (sum_j a_ij Psi_j(u)),
and gamma(t) = |t| * prod over subtrees gamma(u).
"""

from __future__ import annotations

from fractions import Fraction
from itertools import combinations_with_replacement
from typing import Sequence

Tree = tuple  # a root plus a tuple of subtrees; () is the single node


def trees_of_order(n: int) -> list[Tree]:
    """All rooted trees with exactly n nodes, up to isomorphism."""
    if n == 1:
        return [()]
    out = []
    # Partition the n-1 non-root nodes among an unordered multiset of subtrees.
    def partitions(total, max_order):
        if total == 0:
            yield ()
            return
        for k in range(min(total, max_order), 0, -1):
            for rest in partitions(total - k, k):
                yield (k,) + rest
    for part in partitions(n - 1, n - 1):
        groups = []
        for order in sorted(set(part)):
            mult = part.count(order)
            groups.append(list(combinations_with_replacement(trees_of_order(order), mult)))
        stack = [()]
        for grp in groups:
            stack = [acc + choice for acc in stack for choice in grp]
        out.extend(tuple(sorted(s, key=repr)) for s in stack)
    return sorted(set(out), key=repr)


def gamma(t: Tree) -> int:
    """The density gamma(t): |t| times the product of the subtree densities."""
    return order(t) * _prod(gamma(u) for u in t)


def order(t: Tree) -> int:
    return 1 + sum(order(u) for u in t)


def _prod(it):
    r = 1
    for v in it:
        r *= v
    return r


def psi(t: Tree, A: Sequence[Sequence[Fraction]]) -> list[Fraction]:
    """Psi_i(t) for every stage i."""
    s = len(A)
    if not t:
        return [Fraction(1)] * s
    result = [Fraction(1)] * s
    for u in t:
        pu = psi(u, A)
        result = [result[i] * sum(A[i][j] * pu[j] for j in range(s)) for i in range(s)]
    return result


def _square(A):
    """Accept a strictly-lower-triangular tableau in either padded or ragged form."""
    s = len(A)
    return [[Fraction(row[j]) if j < len(row) else Fraction(0) for j in range(s)] for row in A]


def order_conditions(A, b, max_order: int):
    """Yield (tree, order, lhs, rhs, ok) for every condition up to `max_order`."""
    A = _square(A)
    b = [Fraction(x) for x in b]
    for n in range(1, max_order + 1):
        for t in trees_of_order(n):
            lhs = sum(b[i] * p for i, p in enumerate(psi(t, A)))
            rhs = Fraction(1, gamma(t))
            yield t, n, lhs, rhs, lhs == rhs


def effective_order(A, b, max_check: int = 7) -> int:
    """The largest p for which every order condition up to p holds exactly."""
    p = 0
    for n in range(1, max_check + 1):
        if all(ok for _, k, _, _, ok in order_conditions(A, b, n) if k == n):
            p = n
        else:
            break
    return p


def stage_consistency(A, c=None):
    """Check sum_m a[s][m] == c[s] for every stage.

    Violating this makes the stage evaluate the right-hand side at a point that
    is not the one the tableau's c value claims, which is exactly how the
    46732/5147 typo destroyed the method's order.  It is checked separately
    because it localises the fault to a single row.
    """
    A = _square(A)
    out = []
    for s, row in enumerate(A):
        got = sum(row)
        want = Fraction(c[s]) if c is not None else got
        out.append((s, got, want, got == want))
    return out
