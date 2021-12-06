from itertools import product


def get_polar_form(Q):
    return Q + Q.transpose()


class UOV():
    """A class for the (Unbalanced) Oil and Vinegar scheme."""

    def __init__(self, q, m, n, debug=True):
        assert m < n
        self.debug = debug
        self.q = q
        F = GF(q)
        self.F = F
        self.m = m
        self.n = n
        self.V = VectorSpace(F, n)
        self.W = VectorSpace(F, m)
        self.R = PolynomialRing(F, ['x%s' % p for p in range(
            1, n + 1)])
        self.R.inject_variables()
        self.xx = vector(self.R.gens()[:n])
        self.yy = vector(self.R.gens()[n:])
        self.FF = self.construct_central_map()
        self.T, self.PP, self.MM = self.hide_central_map()
        self.O = self.find_oil_subspace()

    def construct_central_map(self):
        m, n = self.m, self.n
        FF = []
        for _ in range(m):
            Q = Matrix(self.F, n)
            for i in range(n):
                for j in range(n):
                    if i >= n - m and j >= n - m:
                        continue
                    Q[i, j] = self.F.random_element()
            FF.append(Q)
        return FF

    def hide_central_map(self):
        m, n = self.m, self.n
        while True:
            T = random_matrix(self.F, n)
            if T.rank() == n:
                break
        PP = [T.transpose() * Q * T for Q in self.FF]
        MM = [get_polar_form(P) for P in PP]
        assert MM == [get_polar_form(P) for P in PP]
        return T, PP, MM

    def find_oil_subspace(self):
        m, n = self.m, self.n
        O2_basis = [self.V([0] * n) for _ in range(m)]
        for i in range(m):
            O2_basis[i][i + n - m] = 1
        O_basis = [self.T.inverse() * o for o in O2_basis]
        O = self.V.subspace(O_basis)

        if self.debug:
            for o in O2_basis:
                for Q in self.FF:
                    assert o * Q * o == 0
            for o in O_basis:
                for P in self.PP:
                    assert o * P * o == 0
            for i in range(m):
                for j in range(i, m):
                    for M in self.MM:
                        assert O_basis[i] * M * O_basis[j] == 0
        return O

    def kipnis_shamir(self, verbose=False):
        m, n = self.m, self.n
        assert n == 2 * m
        inv_subspaces = []
        for i in range(m):
            for j in range(i + 1, m):
                try:
                    Mij = self.MM[i] * self.MM[j].inverse()
                except ZeroDivisionError:
                    continue
                if verbose:
                    print("i,j:", i, j)
                    print("Mi,Mj:", self.MM[i], self.MM[j])
                    print("Mij:", Mij)
                poly = Mij.characteristic_polynomial().radical()
                if poly(Mij) == 0:
                    continue
                inv_subspaces.append(poly(Mij).kernel())
        return inv_subspaces

    def intersection_attack(self, verbose=False):
        m, n = self.m, self.n
        assert n < 3 * m
        found = 0
        equations = []
        for i in range(m):
            for j in range(i + 1, m):
                try:
                    Mi_inv = self.MM[i].inverse()
                    Mj_inv = self.MM[j].inverse()
                    found = 1
                except ZeroDivisionError:
                    continue
                if verbose:
                    print("i,j:", i, j)
                    print("Mi_inv, Mj_inv:\n", Mi_inv,
                          "\n\n", Mj_inv, Mj_inv * self.xx)
                if found == 1:
                    u = Mi_inv * self.xx
                    v = Mj_inv * self.xx
                    for P in self.PP:
                        equations.append(u * P * u)
                        equations.append(v * P * v)
                    for M in self.MM:
                        equations.append(u * M * v)
                    return (equations, i, j)
                return ([], i, j)

    def intersection_attack_advanced(self, verbose=False):
        m, n = self.m, self.n
        assert n > 2 * m and n < 3 * m
        k = 2
        while True:
            if verbose:
                print(k, RR(n / m), RR((2 * k - 1) / (k - 1)))
            if n >= (2 * k - 1) / (k - 1) * m:
                k -= 1
                break
            k += 1
        if verbose:
            print("k:", k)
        LL = []
        for i in range(k):
            L = Matrix(self.F, n)
            while not L.is_invertible():
                for j in range(len(self.MM)):
                    L += self.F.random_element() * self.MM[j]
            LL.append(L)
        if verbose:
            print(LL)
        equations = []
        for i in range(k):
            u = LL[i].inverse() * self.xx
            for P in self.PP:
                equations.append(u * P * u)
            for j in range(i + 1, k):
                v = LL[j].inverse() * self.xx
                for M in self.MM:
                    equations.append(u * M * v)
        return (equations, None, None)


def guess_solve(equations, m, n, xx):
    assert len(equations) == 3 * m
    head = 2 * n - 3 * m
    tail = 3 * m - n
    solution = [1] * tail
    for eq in equations:
        eq = eq(*xx[:head], *([1] * tail))
    for guesses in product(*([GF(q)] * head)):
        solved = 0
        for eq in equations:
            if eq(*guesses, *([1] * tail)) == 0:
                solved += 1
            else:
                continue
        if solved == len(equations):
            return vector(list(guesses) + solution)
    return vector([])

# def guess_solve_advanced(equations, m, n, k, xx):
#     assert len(equations) == m * k * (k + 1) / 2
#     head = k * n - (2 * k - 1) * m
#     tail = k * m - (k - 1) * (n - m)
#     print(len(equations) - head - tail)
#     solution = [1] * tail
#     for eq in equations:
#         eq = eq(*xx[:head], *([1] * tail), *xx)
#     for guesses in product(*([GF(q)] * head)):
#         solved = 0
#         for eq in equations:
#             if eq(*guesses, *([1] * tail), *xx) == 0:
#                 solved += 1
#             else:
#                 continue
#         if solved == len(equations):
#             return vector(list(guesses) + solution)
#     return vector([])

def check_solution(equations, solution):
    for eq in equations:
        assert eq(*solution) == 0
    print("The solution is correct")


q = 11
m = 3
n = 6
uov = UOV(q, m, n)
xx = uov.xx
equations, i, j = uov.intersection_attack()
print("The system to be solved:")
for eq in equations:
    print(eq)
solution = guess_solve(equations, m, n, xx)
if solution == vector([]):
    print("No solution found")
else:
    print(solution)
    check_solution(equations, solution)
    if uov.MM[i].inverse() * uov.V(solution) in uov.O:
        print("The solution corresponds to a vector in O")

    if uov.MM[j].inverse() * uov.V(solution) in uov.O:
        print("The solution corresponds to a vector in O")
