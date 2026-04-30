from power_residues import scan, verify, represent
from power_residues.periods import _constant_term, _is_eth_power

def test_scan_quadratic_small():
    # Quadratic residues mod small primes — classically known
    result = scan(2, 50)
    assert 7 in result   # 2 is a QR mod 7
    assert 17 in result  # 2 is a QR mod 17
    assert 3 not in result  # 3 ≡ 1 mod 2 fails the condition

def test_constant_term_known():
    # q = 7, e = 2: periods are (ζ² + ζ⁴ + ζ⁶) and (ζ + ζ³ + ζ⁵)
    # minimal poly is x² + x − 2, constant term is −2
    assert _constant_term(7, 2) == -2

def test_represent_sum_of_squares():
    # Fermat: p ≡ 1 mod 4 iff p = a² + b²
    eq = lambda a, b: a**2 + b**2
    assert represent(5, eq) == (1, 2)
    assert represent(13, eq) == (2, 3)
    assert represent(3, eq) is None
