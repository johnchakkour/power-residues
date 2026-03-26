"""
This is a Python implementation of calculations from the paper

    Thaine, F (1996) "Properties that characterize Gaussian periods
    and cyclotomic numbers" Proc. Amer. Math. Soc. 124 (1): 35–45.

Let q = ef + 1 be prime, ζ_q a primitive q-th root of unity, and
η_0, ..., η_{e–1} the periods of degree e of Q(ζ_q). The periods
are conjugate and have the same minimal polynomial, which we can
obtain as the characteristic polynomial of the matrix A=[a{i,j}]
where η_0 * η_i = Σ^{e–1}_{j=0} (a{i,j} * η_j). This irreducible
polynomial of degree e is equal to x^e – tr(η_0) + ··· + N(η_0),
where tr and N are respectively the field trace and norm of η_0.
Moreover, note that tr(η_0) = tr(A) = Σ_i η_i = –1, and N(η_0) =
det(A) = (–1)^e * Π_i η_i, so the characteristic polynomial of A
is equal to x^e + x^{e–1} + ··· + (–1)^e * Π_i η_i.
"""

from sympy import *
import cmath


def subgroup(q, prim, e, f):
    # prim is a primitive root mod q
    # c0 is the subgroup of e-th power residues mod q, of order f
    c0 = [pow(prim, e*k, q) for k in range(f)]
    # c1, c2, ..., c_{e-1} are its cosets
    cosets = [[(pow(prim, k, q) * i) % q for i in c0] for k in range(e)]
    return cosets


def periods(sym, sub):
    # Returns the Gaussian period of the subgroup sub, symbolically
    return sum(sym**i for i in sub)


def period_value(sub, q):
    # Returns the Gaussian period of the subgroup sub, algebraically
    z = cmath.exp(2j * cmath.pi / q)
    return sum(z**h for h in sub)


def reduce_powers(expr, z, n):
    # Reduces the exponents in the Gaussian period mod n
    expr = expr.replace(
        lambda x: x.is_Pow and x.base == z,
        lambda x: z**(x.exp % n)
    )
    expr = expand(expr)
    # Use 1 + ζ + ζ^2 + ... + ζ^(n-1) = 0 to eliminate constant term
    const = expr.as_coefficients_dict().get(1, 0)
    if const:
        expr = expand(expr - const * (1 + sum(z**k for k in range(1, n))))
    return expr


def linear(basis_expr, expr, sym):
    # Pick one exponent from each period and read off its coefficient
    def first_exp(e):
        for term in Add.make_args(e):
            if term.is_Pow:
                return term.exp
            if term == sym:
                return 1
    return [expr.coeff(sym, first_exp(b)) for b in basis_expr]
