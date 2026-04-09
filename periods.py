from sympy import *
from typing import *


def _subgroup(q: int, prim: int, e: int) -> list[list[int]]:
    """Returns the subgroup of e-th power residues mod q and its cosets"""
    f = (q - 1) // e
    # prim is a primitive root mod q
    # c0 is the subgroup of e-th power residues mod q, of order f
    c0 = [pow(prim, e * k, q) for k in range(f)]
    # c1, c2, ..., c_{e-1} are its cosets
    cosets = [[(pow(prim, k, q) * i) % q for i in c0] for k in range(e)]
    return cosets


def _periods(sym, sub: list[int]):
    """Returns the Gaussian period of the subgroup sub, symbolically"""
    # Example: if sym = Symbol('ζ') and sub = [2, 4, 6], then
    # _periods(sym, sub) returns ζ**6 + ζ**4 + ζ**2
    return sum(sym ** i for i in sub)


def _reduce_powers(expr, z, n: int):
    """Reduces the exponents in the Gaussian period mod n"""
    expr = expr.replace(
        lambda x: x.is_Pow and x.base == z,
        lambda x: z ** (x.exp % n))
    expr = expand(expr)
    # Use 1 + ζ + ζ^2 + ... + ζ^(n-1) = 0 to eliminate constant term
    const = expr.as_coefficients_dict().get(1, 0)
    if const:
        expr = expand(expr - const * (1 + sum(z ** k for k in range(1, n))))
    return expr


def _linear(basis_expr, expr, sym):
    """Pick one exponent from each period and read off its coefficient"""
    def first_exp(e):
        for term in Add.make_args(e):
            if term.is_Pow:
                return term.exp
            if term == sym:
                return 1
    return [expr.coeff(sym, first_exp(b)) for b in basis_expr]


def _constant_term(q: int, e: int) -> int:
    """Returns the constant term of the minimal polynomial of the Gaussian
     period for q = ef + 1
    """
    g = primitive_root(q)
    cosets = _subgroup(q, g, e)
    z = symbols('ζ')
    periods_list = [_periods(z, ck) for ck in cosets]
    rows = []
    for i in range(e):
        prd = _reduce_powers(
            expand(periods_list[0] * periods_list[i]), z, q)
        rows.append(_linear(periods_list, prd, z))
    mat = Matrix(rows)
    lam = symbols('λ')
    # Constant term is the value of the characteristic polynomial at λ = 0
    return int(mat.charpoly(lam).as_expr().subs(lam, 0))


def _is_eth_power(n: int, q: int, e: int) -> bool:
    """Return True if n is an e-th power residue mod q"""
    return pow(n, (q - 1) // e, q) == 1


def scan(e: int, bound: int, prnt=False) -> list[int]:
    """Scan all primes q = ef + 1 up to bound and print a comparison table."""
    if prnt:
        print(f"e = {e}, scanning primes q = e·f + 1 up to {bound}")
        print()
        print(f"{'q':>6}  {'f':>6}  {'const term':>12}  {'even?':>6}  "
              f"{'2 e-th power?':>14}  {'agree?':>7}")
        print("-" * 62)
    q = 2
    total = 0  # number of primes q such that q – 1 = 0 (mod e)
    res_lst = []  # list of primes q for which 2 is an e-th power mod q
    while q <= bound:
        if (q - 1) % e == 0:
            total += 1
            ct = _constant_term(q, e)
            even = ct % 2 == 0
            if even:
                res_lst.append(q)
            eth_power = _is_eth_power(2, q, e)
            agree = "✓" if even == eth_power else "✗ MISMATCH"
            f = (q - 1) // e
            if prnt:
                print(f"{q:>6}  {f:>6}  {ct:>12}  "
                      f"{'Yes' if even else ' –':>6} {str(eth_power):>14}"
                      f"  {agree:>7}")
        q = nextprime(q)
    if prnt:
        print()
        if e == 2:
            print(f"The primes q <= {bound} for which 2 is a quadratic "
                  f"residue mod q are: ")
        elif e == 3:
            print(f"The primes q <= {bound} for which 2 is a cubic residue "
                  f"mod q are: ")
        else:
            print(f"The primes q <= {bound} for which 2 is a {e}-th power "
                  f"residue mod q are: ")
        for p in res_lst:
            print("    ", end="")
            print(p)
        print()
        print(f"The proportion of primes q <= {bound} satisfying this "
              f"condition is {round(len(res_lst)/total, 3)}.")
    return res_lst


def represent(q: int, eq: Callable[..., int]) -> tuple[int | None, ...] | None:
    """Returns the first tuple of non-negative integers satisfying
    eq(*args) == q, or None if none is found. Assumes eq is non-decreasing
    in each argument.
    """
    a = 0
    while eq(a, 0) <= q:
        b = 0
        while eq(a, b) <= q:
            if eq(a, b) == q:
                return a, b
            b += 1
        a += 1
    return None


def verify(n: int, eq: Callable[..., int], target: int,
           eq_label: str = "q = f(a,b)") -> None:
    """Verify that 2 is an n-th power residue mod q iff eq(a, b) == q has
    a solution, for the first `target` primes q ≡ 1 (mod n).
    """
    header = f"2^((q-1)/{n})≡1"
    print(f"Checking: {eq_label}")
    print(f"{'q':>10}  {header:>15}  {'a':>6}  {'b':>6}")
    print("-" * 46)
    q, count = 2, 0
    failures = []
    while count < target:
        q = nextprime(q)
        if q % n != 1:
            continue
        count += 1
        res = _is_eth_power(2, q, n)
        res_str = "True" if res else "-"
        sol = represent(q, eq) if res else None
        a_str = str(sol[0]) if sol else "-"
        b_str = str(sol[1]) if sol else "-"
        print(f"{q:>10}  {res_str:>15}  {a_str:>6}  {b_str:>6}")
        if res != (sol is not None):
            failures.append(q)
    print("-" * 46)
    if failures:
        print(f"EQUIVALENCE FAILED for: {failures}")
    else:
        print(f"Equivalence verified for all {target} primes q ≡ 1 (mod {n}).")
