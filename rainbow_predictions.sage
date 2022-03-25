'''Adapted from https://github.com/WardBeullens/BreakingRainbow'''

import pandas as pd


def main(q, n, m, o2):
    M = m
    N = n - m
    if q % 2 == 0:
        M -= 1
        N -= 1
    R = PowerSeriesRing(ZZ, 'X')
    X = R.gen()

    def ps_monomials(M, N):
        return 1 / (1 - X) ^ N

    def ps_reg(M, N):
        return (1 - X ^ 2) ^ M / (1 - X) ^ N

    def ps_ranks(M, N):
        return (1 - (1 - X ^ 2) ^ M) / (1 - X) ^ N

    def delete_powers(eq):
        return sum([radical(mon) for mon in eq.monomials()])

    print("M = {}, N = {}\nMacaulay matrices at degree D:".format(M, N))
    # find the first non-positive coefficient in the Hilbert series
    L = 2
    while True:
        L += 1
        if ps_reg(M, N).coefficients()[L - 1] <= 0:
            break

    NumberOfMonomials = ps_monomials(M, N).coefficients()[:L]
    Expected_Coranks = ps_reg(M, N).coefficients()[:L]
    Expected_Ranks = [n - c for n,
                      c in zip(NumberOfMonomials, Expected_Coranks)]
    NumberOfRows = [M * binomial(N + D - 3, N - 1) for D in range(L)]
    data = [NumberOfMonomials[2:], Expected_Ranks[2:], NumberOfRows[2:]]
    df = pd.DataFrame(data, columns=range(2, L))
    df.index = ["Number of cols/mons:", "Expected ranks:", "Number of rows:"]
    print(df.to_string())

    def multiplications(D, N):
        return 3 * binomial(N - 1 + D, D) ^ 2 * binomial(N + 1, 2)

    def complexity(D, N, M, q):
        # we assume that the system is homogenous, so the first guess is free,
        # but we do not account for this in either variant
        guesses = find_number_of_guesses(D, N, M)
        without_guesses = multiplications(D, N).nbits()
        with_guesses = multiplications(
            D, N - guesses).nbits() + ceil(guesses * log(q, 2))
        return (without_guesses, with_guesses)

    def find_number_of_guesses(D, N, M):
        assert D >= 2
        exp_rank = Expected_Ranks[D]
        number_of_guesses = 0
        while exp_rank < binomial(N - number_of_guesses + D, D):
            number_of_guesses += 1
        return number_of_guesses

    print("\nCost of linearization at degree D (in bits):")
    complexities = [complexity(D, N, M, q) for D in range(2, L)]
    data = list(zip(*complexities))
    df = pd.DataFrame(data, columns=range(2, L))
    df.index = ["Without guessing:   ", "With guessing:"]
    print(df.to_string())


# Table 2
# q, n, m, o2 = 31, 30, 20, 10
# q, n, m, o2 = 31, 45, 30, 15
# q, n, m, o2 = 31, 60, 40, 20

# Table 3
# q, n, m, o2 = 16, 30, 20, 10
# q, n, m, o2 = 16, 36, 24, 12
# q, n, m, o2 = 16, 42, 28, 14

# q, n, m, o2 = 2, 36, 24, 12
# q, n, m, o2 = 2, 42, 28, 14
# q, n, m, o2 = 2, 15, 10, 5

# NIST SL1
q, n, m, o2 = 16, 96, 64, 32
# q, n, m, o2 = 2, 96, 64, 32

if __name__ == '__main__':
    main(q, n, m, o2)
