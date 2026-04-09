from sympy import *
from periods import _subgroup, _periods, _linear, _reduce_powers

if __name__ == '__main__':
    try:
        q = int(input("Enter a prime number q: "))
        while not isprime(q):
            q = int(input(f"q = {q} is not prime. Try again: "))
        print(f"q = {q} is prime")
        print()

        print(f"Here are the prime factors of q – 1:")
        factors = list(factorint(q-1).keys())
        print("    ", end="")
        print(*factors, sep=", ")
        e = int(input("Choose a prime e from this list: "))
        while e not in factors:
            e = int(input(f"{e} is not a prime factor of q – 1. Try again: "))
        print()

        print(f"We have q = ef + 1 = {e}·{(q-1) // e} + 1")
        print()

        g = primitive_root(q)
        print(f"g = {g} is a primitive root mod {q}")
        print()

        cosets = _subgroup(q, g, e)
        print(f"In (Z/{q}Z)* we have the subgroup C_0 of e-th power residues"
              f" mod q, and its cosets C_k = g^k·C_0:")
        for k, ck in enumerate(cosets):
            print(f"    C_{k} = {ck}")
        print()

        zeta = Symbol('ζ')
        periods_list = [_periods(zeta, ck) for ck in cosets]
        eta = symbols(f'η0:{e}')  # creates η0, η1, ..., η(e-1)
        print(f"From these we get the following Gaussian periods:")
        for k, pk in enumerate(periods_list):
            print(f"    η{k} = {pk}")
        print()

        print("Products of η0 with each period:")
        rows = []
        for i in range(e):
            prod = _reduce_powers(
                expand(periods_list[0] * periods_list[i]), zeta, q)
            coeffs = _linear(periods_list, prod, zeta)
            rows.append(coeffs)
            lhs = f"η0² " if i == 0 else f"η0·η{i}"
            combo = " + ".join(f"({c})·η{j}" for j, c in enumerate(coeffs))
            print(f"    {lhs} = {combo}")
        print()

        lam = symbols('λ')
        A = Matrix(rows)
        print("The matrix of linear coefficients is")
        print()
        pprint(A)
        print()
        print("with characteristic polynomial")
        print()
        char_poly = A.charpoly(lam).as_expr()
        pprint(factor(char_poly))
        print()
        if e == 2:
            print("η0, η1", end="")
        elif e == 3:
            print("η0, η1, η2", end="")
        else:
            print(f"η0, ..., η{e-1}", end="")
        print(f" generate the unique subfield of degree {e} over Q in"
              f" Q({zeta}_{q}).")

    except ValueError:
        print("Please enter a valid integer.")
