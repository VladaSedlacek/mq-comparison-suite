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
# q, n, m, o2 = 16, 36, 24, 12
# q, n, m, o2 = 16, 42, 28, 14

# q, n, m, o2 = 2, 36, 24, 12
# q, n, m, o2 = 2, 42, 28, 14
# q, n, m, o2 = 2, 15, 10, 5

q, n, m, o2 = 2, 12, 8, 4
# q, n, m, o2 = 31, 18, 12, 6
# q, n, m, o2 = 4, 9, 6, 3

# NIST SL1
# q, n, m, o2 = 16, 96, 64, 32

seed = 0
set_random_seed(seed)

K = GF(q)
z = K.gens()[0]
ext_deg = K.modulus().degree()
# R = K
R = PolynomialRing(K, 'X')
X = R.gen()

attempts = 0

basis_Fn = (R**n).basis()
basis_Fm = (R**m).basis()

# if you run without O2, the system is not guaranteed to have a solution


def Attack(PK, O2=None):
    global attempts

    # pick a random vector x
    x = vector([K.random_element() for i in range(n)])
    x = vector([R(e) for e in x])
    # x[-1] = X
    print("x:", x)
    while Eval(PK, x)[0] == 0:
        x = vector([K.random_element() for i in range(n)])

    # compute linear map D_x = P'(x,.)
    D_x = Matrix(R, [Differential(PK, x, b) for b in basis_Fn])
    # print(D_x)
    D_x_ker = Matrix(D_x.kernel().basis())
    print("REF:")
    print(D_x_ker())
    print("RREF:")
    print(D_x_ker.rref())

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
SS = Attack(PK)
# print(SS[0])
# SS = [ Matrix(N,N,[ K.random_element() for _ in range(N*N) ])  for _ in SS]

N = SS[0].ncols()
M = len(SS)
PR = PolynomialRing(R, N, 'x')
PR.inject_variables(verbose=False)
x_vec = vector(PR, [PR.gens()])
equations_orig = [x_vec * M * x_vec for M in SS]


def linear_combination(coefficients, objects):
    assert len(coefficients) == len(objects)
    return sum(c * o for c, o in zip(coefficients, objects))


def weil_decomposition(poly):
    if poly == 0:
        return []
    # Constant coefficients come first
    extension_coeffs = [c.polynomial().list() for c in poly.coefficients()]
    max_len = max(len(ec) for ec in extension_coeffs)
    # Pad the coefficient lists
    for ec in extension_coeffs:
        for _ in range(max_len - len(ec)):
            ec.append(0)
    base_coeff_list = zip(*extension_coeffs)
    base_polys = [linear_combination(coeffs, poly.monomials())
                  for coeffs in base_coeff_list]
    return base_polys


weil_ring = PolynomialRing(
    K, ['w%s_%s' % (p1, p2) for p1, p2 in itertools.product(range(1, n - m + 1), range(ext_deg))], order='degrevlex')
weil_ring.inject_variables(verbose=False)
ww = weil_ring.gens()
ww_parts = [ww[i:i + ext_deg] for i in range(0, ext_deg * (n - m), ext_deg)]
zs = [z ^ i for i in range(ext_deg)]
yy_weil = vector([linear_combination(w, zs)
                  for w in ww_parts])
yy_weil_affine = vector(list(yy_weil[1:-1]) + [1])
equations = [yy_weil_affine * s * yy_weil_affine for s in SS]
equations_weil = [
    w_eq for eq in equations for w_eq in weil_decomposition(eq)]


print("Original:")
for eq in equations_orig[:1]:
    print(eq)

# print("Weil:")
# for eq in equations_weil[:ext_deg]:
#     print(eq)

# print(K.modulus())


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


# for D in range(2, 10):
#     eqns = []
#     for p in equations:
#         for Mon in Monomials(PR.gens(), D - 2):
#             NewMon = (p * Mon)
#             NewMon = NewMon(x0=0, x1=1)
#             if NewMon != 0:
#                 if q == 2:
#                     NewMon = delete_powers(NewMon)
#                 eqns.append(NewMon)
#     s = Sequence(eqns)
#     Matrix, Mon = s.coefficient_matrix()
#     linear = [PR(m) for m in Mon if PR(m).degree() <= 1]
#     rank = Matrix.rank()
#     rows = Matrix.nrows()
#     cols = Matrix.ncols()
#     print("\nD = ", D)
#     print("\t\trank: %d,\t cols: %d, rows: %d" % (rank, cols, rows))
#     print("\tbefore guessing:\t cols: %d, rows: %d" %
#           (binomial(N - 1 + D, D), M * binomial(N + D - 3, N - 1)))
#     if rank == cols - 1:
#         print("SOLVABLE!")
#         linear = [PR(m) for m in Mon if PR(m).degree() <= 1]
#         print(linear)
#         for bv in kernel(Matrix.transpose()).basis():
#             print(bv[-len(linear):])
