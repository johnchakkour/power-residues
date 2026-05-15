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


def _periods(sym, sub: list[int]) -> Expr:
    """Returns the Gaussian period of the subgroup sub, symbolically"""
    # Example: if sym = Symbol('ζ') and sub = [2, 4, 6], then
    # _periods(sym, sub) returns ζ**6 + ζ**4 + ζ**2
    return Add(*[sym ** i for i in sub])


def _reduce_powers(expr: Expr, z, n: int) -> Expr:
    """Reduces the exponents in the Gaussian period mod n"""
    d = expr.as_coefficients_dict()
    red: dict[int, int] = {}
    for term, coef in d.items():
        expr = 0 if term == 1 else (term.exp if term.is_Pow else 1)
        red[expr % n] = red.get(expr % n, 0) + int(coef)
    # Eliminate constant via 1 + ζ + ... + ζ^(n-1) = 0
    const = red.pop(0, 0)
    if const:
        for k in range(1, n):
            red[k] = red.get(k, 0) - const
    return Add(*[c * z**e for e, c in red.items() if c])


def _linear(basis_expr, expr, sym) -> list[int]:
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
    # lam = symbols('λ')
    # Constant term is the value of the characteristic polynomial at λ = 0
    # return int(mat.charpoly(lam).as_expr().subs(lam, 0))
    return int(((-1) ** e) * mat.det())


def _is_eth_power(n: int, q: int, e: int) -> bool:
    """Return True if n is an e-th power residue mod q"""
    return pow(n, (q - 1) // e, q) == 1


def _prompt_prime(prompt: str, reprompt: str) -> int:
    """Prompt the user until they enter a prime"""
    n = int(input(prompt))
    while not isprime(n):
        n = int(input(reprompt.format(n=n)))
    return n


def _prompt_choice(prompt: str, reprompt: str, choices: list) -> int:
    """Prompt the user until they enter a value from choices"""
    v = int(input(prompt))
    while v not in choices:
        v = int(input(reprompt.format(v=v)))
    return v


def _first_row(q: int, e: int) -> list[int]:
    """Coefficients (a_{0,0}, ..., a_{0,p-1}) of η0² = Σ_i a_{0,i}·η_i."""
    g = primitive_root(q)
    cosets = _subgroup(q, g, e)
    z = symbols('ζ')
    periods_list = [_periods(z, ck) for ck in cosets]
    prd = _reduce_powers(expand(periods_list[0] ** 2), z, q)
    return [int(c) for c in _linear(periods_list, prd, z)]


def scan(e: int, bound: int, prnt=False) -> list[int]:
    """Scan all primes q = ef + 1 up to bound and print a comparison table."""
    if prnt:
        print(f"We have e = {e}. For which primes q is 2 an e-th power "
              f"residue mod q?")
        print(f"Note that necessarily such primes must satisfy "
              f"q ≡ 1 (mod {e}).")
        print()
        print("We suspect that an equivalent condition to 2 being an e-th "
              "power residue")
        print("mod q is that the constant term of the minimal polynomial of "
              "the Gaussian")
        print("period η0 associated to the subgroup of eth power residues "
              "in (Z/qZ)* is even.")
        print()
        print(f"Scanning primes q = {e}·f + 1 up to {bound}:")
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
        elif e == 5:
            print(f"The primes q <= {bound} for which 2 is a quintic residue "
                  f"mod q are: ")
        else:
            print(f"The primes q <= {bound} for which 2 is a {e}-th power "
                  f"residue mod q are: ")
        for p in res_lst:
            print("    ", end="")
            print(p)
        print()
        if total:
            print(f"The proportion of primes q ≡ 1 (mod {e}) <= {bound} "
                  f"satisfying this condition is "
                  f"{round(len(res_lst)/total, 3)}.")
        print(f"The expected proportion of primes satisfying this condition "
              f"is 1/e = 1/{e} = {round(1/e, 3)}.")
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
        sol = represent(q, eq)
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


def gaussian_periods(q: int | None = None, e: int | None = None) -> None:
    """
    For a prime q = ef + 1, display the subgroup C of e-th power residues in
    (Z/qZ)*, the Gaussian periods associated to C and its e – 1 cosets, their
    pairwise products expressed as linear combinations of the periods, and
    their resulting minimal polynomial over Q.

    Parameters
    ----------
    q : prime of the form ef + 1; prompted interactively if None.
    e : prime factor of q – 1 determining the subgroup index; prompted if None.
    """
    try:
        if q is None:
            q = _prompt_prime(
                "Enter a prime number q: ",
                "q = {n} is not prime. Try again: ")
        elif not isprime(q):
            raise ValueError(f"{q} is not prime.")

        print(f"q = {q} is prime\n")

        factors = list(factorint(q - 1).keys())
        print("Here are the prime factors of q – 1:")
        print("    " + ", ".join(str(p) for p in factors) + "\n")

        if e is None:
            e = _prompt_choice(
                "Choose a prime e from this list: ",
                "{v} is not a prime factor of q – 1. Try again: ",
                factors)
        elif e not in factors:
            raise ValueError(f"{e} is not a prime factor of q – 1 = {q - 1}.")

        f = (q - 1) // e
        print(f"We have q = ef + 1 = {e}·{f} + 1\n")

        g = primitive_root(q)
        print(f"g = {g} is a primitive root mod {q}\n")

        cosets = _subgroup(q, g, e)
        print(f"In (Z/{q}Z)* we have the subgroup C_0 of e-th power residues "
              f"mod q, and its cosets C_k = g^k·C_0:")
        for k, ck in enumerate(cosets):
            print(f"    C_{k} = {ck}")
        print()

        zta = Symbol('ζ')
        periods_list = [_periods(zta, ck) for ck in cosets]
        print("From these we get the following Gaussian periods:")
        for k, pk in enumerate(periods_list):
            print(f"    η{k} = {pk}")
        print()

        print("Products of η0 with each period:")
        rows = []
        for i in range(e):
            prd = _reduce_powers(
                expand(periods_list[0] * periods_list[i]), zta, q)
            coeffs = _linear(periods_list, prd, zta)
            rows.append(coeffs)
            lhs = "η0²" if i == 0 else f"η0·η{i}"
            combo = " + ".join(f"({c})·η{j}" for j, c in enumerate(coeffs))
            print(f"    {lhs} = {combo}")
        print()

        lam = symbols('λ')
        a = Matrix(rows)
        print("The matrix of linear coefficients is\n")
        pprint(a)
        print()

        char_poly = a.charpoly(lam).as_expr()
        print("with characteristic polynomial\n")
        pprint(factor(char_poly))
        print()

        period_str = ", ".join(
            f"η{k}" for k in range(e)) if e <= 3 else f"η0, ..., η{e - 1}"
        print(f"{period_str} generate the unique subfield of degree {e} "
              f"over Q in Q(ζ_{q}).")

    except ValueError as exc:
        print(f"Error: {exc}")


def check_row(e: int, bound: int, prnt: bool = True) -> list[list[int]]:
    """For each prime q = ef + 1 ≤ bound, compute the row
    (a_{0,0}, ..., a_{0,e-1}) such that η0² = Σ_i a_{0,i}·η_i, and bucket
    q by which coefficients are odd.

    Returns a list of e lists: result[i] is the primes q for which a_{0,i}
    is odd.
    """
    if not isprime(e):
        raise ValueError(f"{e} is not prime.")

    odd_primes: list[list[int]] = [[] for _ in range(e)]

    if prnt:
        print(f"Scanning primes q = {e}·f + 1 up to {bound}.")
        print(f"For each q we expand η0² = a_{{0,0}}·η0 + ... "
              f"+ a_{{0,{e-1}}}·η{e-1} and record the parity of each "
              f"coefficient.\n")
        coeff_hdr = "  ".join(f"a_0,{i}".rjust(8) for i in range(e))
        header = f"{'q':>8}  {coeff_hdr}  parity"
        print(header)
        print("-" * len(header))

    q = 2
    while q <= bound:
        if (q - 1) % e == 0:
            row = _first_row(q, e)
            for i, a in enumerate(row):
                if a % 2 != 0:
                    odd_primes[i].append(q)
            if prnt:
                row_str = "  ".join(f"{a:>8}" for a in row)
                parity = " " + "".join("1" if a % 2 else "0" for a in row)
                print(f"{q:>8}  {row_str}  {parity}")
        q = nextprime(q)

    if prnt:
        print()
        for i in range(e):
            print(f"Primes q ≡ 1 (mod {e}), q ≤ {bound}, with a_{{0,{i}}} odd:")
            if odd_primes[i]:
                print("    " + ", ".join(str(qq) for qq in odd_primes[i]))
            else:
                print("    (none)")
            print()

    return odd_primes
