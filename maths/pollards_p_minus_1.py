from __future__ import annotations

from math import gcd, log
from typing import Iterable
from random import randint


def pollard_rho(
    num: int,
    B: int,
    step: int = 1,
    attempts: int = 3,
) -> int | None:
    """
    Use Pollard's p-1 algorithm to return a nontrivial factor of ``num``.
    It's required for ``num`` to be ``B``-power smooth, that is, that
    its prime factors are less than ``B``.
    The returned factor may be composite and require further factorization.
    If the algorithm will return None if it fails to find a factor within
    the specified number of attempts.
    If ``num`` is prime, this algorithm is guaranteed to return None.
    https://en.wikipedia.org/wiki/Pollard%27s_rho_algorithm
    https://en.wikipedia.org/wiki/Pollard%27s_p_%E2%88%92_1_algorithm
    """

    # The algorithm works only for positive integers
    if num < 0:
        raise ValueError("The input value must be a natural number (positive)")

    # There are no primes bellow 2
    if num < 2:
        return None

    # Pollard's Rho algorithm requires a function that returns pseudorandom
    # values between 0 <= X < ``num``.  It doesn't need to be random in the
    # sense that the output value is cryptographically secure or difficult
    # to calculate, it only needs to be random in the sense that all output
    # values should be equally likely to appear.
    # For this reason, Pollard suggested using ``f(x) = (x**2 - 1) % num``
    # However, the success of Pollard's algorithm isn't guaranteed and is
    # determined in part by the initial seed and the chosen random function.
    # To make retries easier, we will instead use ``f(x) = (x**2 + C) % num``
    # where ``C`` is a value that we can modify between each attempt.
    def generate_primes(limit: int) -> Iterable[int]:
        sieve = iter(range(2,limit+1))
        while True:
            try:
                prime = next(sieve)
            except StopIteration:
                break
            yield prime
            not_multiple_of_prime = lambda x, p=prime: x % p != 0
            sieve = filter(not_multiple_of_prime, sieve)

    for _ in range(attempts):
        guess = randint(2,num)
        for prime in generate_primes(B):
            s = int(log(B, prime))
            r = pow(prime, s)
            guess = pow(guess, r, num)

            divisor = gcd(guess - 1, num)
            if divisor != 1 and divisor != num:
                # The divisor is a nontrivial factor of ``num``!
                return divisor

    # We haven't found a divisor within the requested number of attempts.
    # We were unlucky or ``num`` itself is actually prime.
    return None


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "num",
        type=int,
        help="The value to find a divisor of",
    )
    parser.add_argument(
        "B",
        type=int,
        help="The smoothness of the previous number",
        default=200,
    )
    parser.add_argument(
        "--attempts",
        type=int,
        default=3,
        help="The number of attempts before giving up",
    )
    args = parser.parse_args()

    divisor = pollard_rho(args.num, args.B, attempts=args.attempts)
    if divisor is None:
        print(f"{args.num} is probably prime")
    else:
        quotient = args.num // divisor
        print(f"{args.num} = {divisor} * {quotient}")
