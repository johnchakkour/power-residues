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

