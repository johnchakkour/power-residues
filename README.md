# power_residues

A Python library for exploring **power residue criteria** via the theory of
**Gaussian periods** and cyclotomic fields.

For which primes $q$ is $2$ an $e$-th power residue mod $q$? The classical 
answer is that this is controlled by the arithmetic of field extensions of
$\mathbb{Q}$. This library explicitly investigates the case of cyclotomic
fields and makes the connection computational.

## Mathematical background

Let $q=ef+1$ be prime; then the cyclotomic extension $\mathbb{Q}(\zeta_q)/\mathbb{Q}$
has Galois group $(\mathbb{Z}/q\mathbb{Z})^\times$, which has a unique subgroup
$C_0$ of index $e$ and order $f$. To $C_0$ and its cosets $C_1,\dots,C_{e-1}$ we
respectively associate the Gaussian periods
$$\eta_i = \sum_{j\in C_i} \zeta_q^j.$$
These periods $\eta_0,\dots,\eta_{e-1}$ are conjugate algebraic integers, all 
sharing the same minimal polynomial $\mu(x)$ over $\mathbb{Q}$, and they generate 
the unique subfield of $\mathbb{Q}(\zeta_q)$ of degree $e$ over $\mathbb{Q}$.

### Main observation

Numerical experiments suggest that $2$ is an $e$-th power mod $q$ if and only if
the constant term of $\mu(x)$ is even. This library verifies this equivalence
computationally and explores related criteria (such as needing $q=a^2+27b^2$) that
arise in specific cases.

## Installation

Clone the repository and install in editable mode:
```
git clone https://github.com/johnchakkour/power_residues.git
cd power_residues
pip install -e .
```
Requirements: Python 3.11+, [ SymPy ](https://www.sympy.org/en/index.html).

## Usage

### Exploration

`gaussian_periods()` walks through the full period construction for a prime
$q$ of your choice, displaying the cosets, periods, their product table, and
the resulting minimal polynomial.
```
>>> from power_residues import gaussian_periods
>>> gaussian_periods(q=13, e=3)

q = 13 is prime

Here are the prime factors of q – 1:
    2, 3

We have q = ef + 1 = 3·4 + 1

g = 2 is a primitive root mod 13

In (Z/13Z)* we have the subgroup C_0 of e-th power residues mod q,
and its cosets C_k = g^k·C_0:
    C_0 = [1, 8, 12, 5]
    C_1 = [2, 3, 11, 10]
    C_2 = [4, 6, 9, 7]

...

with characteristic polynomial

 3    2
λ  + λ  - 4·λ + 1

η0, η1, η2 generate the unique subfield of degree 3 over Q in Q(ζ_13).
```
Called without arguments, the function prompts for input interactively.

### Scanning for power residues

`scan(e, bound)` checks all primes $q\equiv 1\pmod{e}$ up to `bound` and compares the
parity of the constant term of $\mu(x)$ against the direct power residue test.
```
from power_residues import scan

# Find all primes q ≤ 200 for which 2 is a cubic residue
result = scan(e=3, bound=200, prnt=True)
```

### Verifying representation criteria

`verify(n, eq, target, eq_label)` checks whether the equivalence "2 is an $n$-th
power residue mod $q\iff q=f(a,b)$ is solvable" holds for the first `target` primes 
$q\equiv 1\pmod{n}$.
```
from power_residues import verify

# 2 is a cubic residue mod q iff q = a² + 27b²
verify(3, lambda a, b: a**2 + 27*b**2, target=20, eq_label="q = a² + 27b²")
```
```
Checking: q = a² + 27b²
         q    2^((q-1)/3)≡1       a       b
----------------------------------------------
         7                -       -       -
        13                -       -       -
        19                -       -       -
        31             True       2       1
        37                -       -       -
        43             True       4       1
        61                -       -       -
----------------------------------------------
Equivalence verified for all 7 primes q ≡ 1 (mod 3).
```

## Project structure
```
power_residues/
├── power_residues/
│   ├── __init__.py       # Public API
│   └── periods.py        # Core implementation
├── tests/
│   └── test_periods.py
├── examples/
│   └── demo.py
├── pyproject.toml
└── README.md
```

## Running tests
```
pip install pytest
pytest
```

## License
MIT License. See [ LICENSE ](LICENSE).
