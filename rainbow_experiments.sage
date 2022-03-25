'''Adapted from https://github.com/WardBeullens/BreakingRainbow'''

import itertools
import pandas as pd

load('rainbow.sage')

# Table 2
# q, n, m, o2 = 31, 30, 20, 10
# q, n, m, o2 = 31, 45, 30, 15
# q, n, m, o2 = 31, 60, 40, 20

# Table 3
# q, n, m, o2 = 16, 30, 20, 10
q, n, m, o2 = 16, 36, 24, 12
# q, n, m, o2 = 16, 42, 28, 14

# q, n, m, o2 = 2, 36, 24, 12
# q, n, m, o2 = 2, 42, 28, 14
# q, n, m, o2 = 2, 15, 10, 5

# NIST SL1
# q, n, m, o2 = 16, 96, 64, 32

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
    if O2 is not None:
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

        print("The expected subsolution is:")
        print(Sol[1:])

    # Compose smaller system D_x(o)= 0 and P(o) = 0
    SS = [D_x_ker * M * D_x_ker.transpose() for M in PK]
    for s in SS:
        Make_UD(s)
    if not O2 is None:
        assert Eval(SS, Sol) == vector([0] * m)

    if q % 2 == 0:
        Px = Eval(PK, x)
        SSS = [(SS[i] * Px[0] + SS[0] * Px[i])[1:, 1:]
               for i in range(1, len(SS))]
        SS = SSS
        if not O2 is None:
            assert Eval(SSS, Sol[1:]) == vector([0] * (m - 1))
    return SS


PK, O2, O1, W = Keygen(q, n, m, o2)
tP = Attack(PK, O2)
# tP = [ Matrix(N,N,[ K.random_element() for _ in range(N*N) ])  for _ in tP]

N = tP[0].ncols()
M = len(tP)
PR = PolynomialRing(K, N, 'x')
PR.inject_variables(verbose=False)
x_vec = vector(PR, [PR.gens()])
tP = [x_vec * M * x_vec for M in tP]


def Monomials(vars, degree):
    '''Compute all monomials of a certain degree'''
    if degree < 0:
        return
    for comb in itertools.combinations_with_replacement(vars, degree):
        u = 1
        for var in comb:
            u *= var
        yield u
    return


R = PowerSeriesRing(ZZ, 'X')
X = R.gen()


def delete_powers(eq):
    return sum([radical(mon) for mon in eq.monomials()])


for D in range(2, 10):
    eqns = []
    for p in tP:
        for Mon in Monomials(PR.gens(), D - 2):
            NewMon = (p * Mon)
            NewMon = NewMon(x0=0, x1=1)
            if NewMon != 0:
                if q == 2:
                    NewMon = delete_powers(NewMon)
                eqns.append(NewMon)
    s = Sequence(eqns)
    Matrix, Mon = s.coefficient_matrix()
    linear = [PR(m) for m in Mon if PR(m).degree() <= 1]
    rank = Matrix.rank()
    rows = Matrix.nrows()
    cols = Matrix.ncols()
    print("\nD = ", D)
    print("\t\trank: %d,\t cols: %d, rows: %d" % (rank, cols, rows))
    print("\tbefore guessing:\t cols: %d, rows: %d" %
          (binomial(N - 1 + D, D), M * binomial(N + D - 3, N - 1)))
    if rank == cols - 1:
        print("SOLVABLE!")
        linear = [PR(m) for m in Mon if PR(m).degree() <= 1]
        print(linear)
        for bv in kernel(Matrix.transpose()).basis():
            print(bv[-len(linear):])
