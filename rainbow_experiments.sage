'''Adapted from https://github.com/WardBeullens/BreakingRainbow'''

import itertools

load('rainbow.sage')

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

seed = 0
set_random_seed(seed)

K = GF(q)
z = K.gens()[0]

attempts = 0

basis_Fn = (K**n).basis()
basis_Fm = (K**m).basis()

# if you run without O2, the system is not guaranteed to have a solution


def Attack(PK, O2=None):
    global attempts

    # pick a random vector x
    x = vector([K.random_element() for i in range(n)])
    while Eval(PK, x)[0] == 0:
        x = vector([K.random_element() for i in range(n)])

    # compute linear map D_x = P'(x,.)
    D_x = Matrix(K, [Differential(PK, x, b) for b in basis_Fn])
    D_x_ker = Matrix(D_x.kernel().basis())

    if q % 2 == 0:
        D_x_ker[0] = x

    if D_x_ker.rank() != n - m:
        return Attack(PK, O2)

    attempts += 1

    Sol = None
    if not O2 is None:
        V = K**n
        I = V.span(D_x_ker).intersection(V.span(O2.transpose()))
        if I.dimension() == 0:
            # print("Attack would fail. resample x")
            return Attack(PK, O2)

        print("Intersection has dimension:", I.dimension())
        Sol = I.basis()[0]

        Sol = D_x_ker.transpose().solve_right(Sol)

        if Sol[-1] == 0:
            print("last entry is zero, resample x")
            return Attack(PK, O2)

        Sol = Sol / Sol[-1]

        print("Good D_x found after %d attempts." % attempts)

        # print("The expected solution is:")
        # print(Sol)
        print("The expected subsolution is:")
        print(Sol[1:])

    # Compose smaller system D_x(o)= 0 and P(o) = 0
    SS = [D_x_ker * M * D_x_ker.transpose() for M in PK]
    for s in SS:
        Make_UD(s)

    assert Sol is not None
    # if not Sol is None:
    #     print("Sanity check: evaluation of sol")
    #     print(Eval(SS, Sol))

    if q % 2 == 0:
        Px = Eval(PK, x)
        SSS = [(SS[i] * Px[0] + SS[0] * Px[i])[1:, 1:]
               for i in range(1, len(SS))]

        # if not Sol is None:
        #     print("Sanity check: evaluation of sol[1:]")
        #     print(Eval(SSS, Sol[1:]))

        SS = SSS

    return SS


PK, O2, O1, W = Keygen(q, n, m, o2)
# print('O2')
# print(O2)
# print('O1')
# print(O1)
# print('W')
# print(W)
tP = Attack(PK, O2)

N = tP[0].ncols()
M = len(tP)

# tP = [ Matrix(N,N,[ K.random_element() for _ in range(N*N) ])  for _ in tP]

print("M = %d, N = %d" % (M, N))
PR = PolynomialRing(K, N, 'x')
PR.inject_variables()

x_vec = vector(PR, [PR.gens()])
tP = [x_vec * M * x_vec for M in tP]

# Compute all monomials of a certain degree


def Monomials(vars, degree):
    if degree < 0:
        return

    for comb in itertools.combinations_with_replacement(vars, degree):
        u = 1
        for var in comb:
            u *= var
        yield u
    return


L = 10
Expected_Ranks = [0] * L
Expected_Ranks[0] = 1

for _ in range(N):
    for i in range(1, L):
        Expected_Ranks[i] += Expected_Ranks[i - 1]

NumberOfMonomials = [x for x in Expected_Ranks]

for _ in range(M):
    for i in range(L - 1, 1, -1):
        Expected_Ranks[i] -= Expected_Ranks[i - 2]

Expected_Ranks = [NumberOfMonomials[i] - Expected_Ranks[i] for i in range(L)]

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


def multiplications(D, N):
    return 3 * binomial(N - 1 + D, D) ^ 2 * binomial(N + 1, 2)


def complexity(D, N, M, q=16):
    # the first guess is free
    guesses = find_number_of_guesses(D, N, M) - 1
    without_guesses = multiplications(D, N).nbits()
    with_guesses = multiplications(
        D, N + 1 - guesses).nbits(), guesses  # * log(q, 2)
    return (without_guesses, with_guesses)


def find_number_of_guesses(D, N, M):
    assert D >= 2
    exp_rank = ps_ranks(M, N).coefficients()[D - 2]
    number_of_guesses = 1
    while exp_rank < binomial(N - number_of_guesses + D, D):
        number_of_guesses += 1
    return number_of_guesses


Expected_Coranks = ps_reg(M, N).coefficients()[:L]
NumberOfRows = [M * binomial(N + D - 3, N - 1) for D in range(L)]

assert NumberOfMonomials == ps_monomials(M, N).coefficients()[:L]
assert Expected_Coranks == ps_reg(M, N).coefficients()[:L]
assert Expected_Ranks == [n - c for n,
                          c in zip(NumberOfMonomials, Expected_Coranks)]
print("Number of columns/monomials:", NumberOfMonomials)
print("Number of rows:             ", NumberOfRows)
# print("Expected co-ranks:          ", Expected_Coranks)
print("Expected ranks (series):    ", Expected_Ranks)
print("Expected ranks (adjusted):  ", [min(rk, r, c) for rk, r, c in zip(
    Expected_Ranks, NumberOfRows, NumberOfMonomials)])


for D in range(2, min(5, N)):
    eqns = []
    for p in tP:
        # print(p(x0=0, x1=0, x2=1, x3=1, x4=0))
        for Mon in Monomials(PR.gens(), D - 2):
            NewMon = (p * Mon)
            # , x1=0, x2=1, x3=0, x4=1, x5=1, x6=0)
            NewMon = NewMon(x0=z + 1, x1=z ^ 3 + z, x2=z ^ 3 +
                            z ^ 2 + 1, x3=z ^ 2 + z + 1, x4=z ^ 3 + z ^ 2 + z, x5=z + 1)

            if NewMon != 0:
                if q == 2:
                    NewMon = delete_powers(NewMon)
                eqns.append(NewMon)

    s = Sequence(eqns)
    M, Mon = s.coefficient_matrix()
    linear = [PR(m) for m in Mon if PR(m).degree() <= 1]
    # print(Mon)
    # print(M.density())
    # rhs = vector(K, [0] * M.nrows())
    # print(parent(M), parent(rhs))
    # print(M.solve_right(rhs))
    Rank = M.rank()
    rows = M.nrows()
    cols = M.ncols()
    print("\nD = ", D)
    print("\trank: %d,\t rows: %d, cols: %d" % (Rank, rows, cols))
    print("\tbefore:\t\t rows: %d, cols: %d" %
          (NumberOfRows[D], NumberOfMonomials[D]))
    if Rank == cols - 1:
        print("SOLVABLE!")
        linear = [PR(m) for m in Mon if PR(m).degree() <= 1]
        print(linear)
        for bv in kernel(M.transpose()).basis():
            print(bv[-len(linear):])
