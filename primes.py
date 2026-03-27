from periods import *

if __name__ == '__main__':
    bound = int(input("Scan primes q up to: "))
    print(f"{'e':>5} "
          f"{'primes q = 1 mod e such that 2 is an e-th power mod q':>90}")
    print("-" * 68)
    for e in primerange(3, 100):
        print(f"{e:>5} {str(scan(e, bound)):>90}")
    """
    try:
        e = int(input("Enter a prime e: "))
        bound = int(input("Scan primes q up to: "))
        print()
        scan(e, bound)
    except ValueError:
        print("Please enter valid integers.")
    """
