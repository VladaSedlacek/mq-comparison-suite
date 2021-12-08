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
        self.reduced = self.q % 2 == 0 and self.n % 2 == 1
        self.V = VectorSpace(F, n)
        self.W = VectorSpace(F, m)
        self.R = PolynomialRing(F, ['x%s' % p for p in range(
            1, n + 1)], order="neglex")
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

    def test_oil_space_membership(self, v):
        return (v in self.O, [v * P * v for P in self.PP])

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

    def intersection_attack(self, verbose=False, advanced=True):
        if self.reduced:
            m, n = self.m - 1, self.n - 1
            MM = [M[1:, 1:] for M in self.MM]
            PP = [P[1:, 1:] for P in self.PP]
            xx = (self.xx)[1:]
        else:
            m, n = self.m, self.n
            MM = self.MM
            PP = self.PP
            xx = self.xx
        if advanced:
            assert n >= 2 * m and n < 3 * m
            k = find_max_k(m, n, True, verbose=verbose)
            LL = []
            for i in range(k):
                while True:
                    coefficients = self.W.random_element()
                    L = linear_combination(coefficients, MM)
                    if L.is_invertible():
                        LL.append(L)
                        break
            equations = []
            for i in range(k):
                u = LL[i].inverse() * xx
                for P in PP:
                    equations.append(u * P * u)
                for j in range(i + 1, k):
                    v = LL[j].inverse() * xx
                    for M in MM:
                        equations.append(u * M * v)
            matrices = [L.inverse() for L in LL]
        else:
            assert n < 3 * m
            found = 0
            equations = []
            for i in range(m):
                for j in range(i + 1, m):
                    try:
                        Mi_inv = MM[i].inverse()
                        Mj_inv = MM[j].inverse()
                        found = 1
                    except ZeroDivisionError:
                        continue
                    if verbose:
                        print("i,j:", i, j)
                        print("Mi_inv, Mj_inv:\n", Mi_inv,
                              "\n\n", Mj_inv, Mj_inv * xx)
                    if found == 1:
                        u = Mi_inv * self.xx
                        v = Mj_inv * self.xx
                        for P in PP:
                            equations.append(u * P * u)
                            equations.append(v * P * v)
                        for M in MM:
                            equations.append(u * M * v)
                        matrices = [Mi_inv, Mj_inv]
        return equations, matrices


def linear_combination(coefficients, objects):
    assert len(coefficients) == len(objects)
    return sum(c * o for c, o in zip(coefficients, objects))


def find_max_k(m, n, verbose=False):
    k = 2
    while True:
        if verbose:
            print("current k:", k, ", n/m:", (n / m).numerical_approx(digits=3),
                  ", (2 * k - 1) / (k - 1):", ((2 * k - 1) / (k - 1)).numerical_approx(digits=3))
        if n >= (2 * k - 1) / (k - 1) * m:
            k -= 1
            break
        if k > max(sqrt(m), 3):
            break
        k += 1
    if verbose:
        print("k:", k)
    return k


def guess_solve(equations, q, m, n, xx, advanced=False, reduced=False, verbose=False):
    if advanced:
        k = find_max_k(m, n, verbose=False)
        total = ZZ(m * k * (k + 1) / 2)
        assert len(equations) == total
        tail = k * m - (k - 1) * (n - m)
    else:
        assert len(equations) == 3 * m
        tail = 3 * m - n
    constraints = [1] * tail
    head = n - tail
    for eq in equations:
        if reduced:
            first_var = 0
        else:
            first_var = xx[0]
        eq = eq(first_var, *xx[1:head], *constraints)
    if reduced:
        guess_space = product(*([GF(q)] * (head - 1)))
    else:
        guess_space = product(*([GF(q)] * head))

    for guesses in guess_space:
        if reduced:
            values = 0, *guesses, *constraints
        else:
            values = *guesses, *constraints
        if verbose:
            print("guess:", values)
        solved = 0
        for eq in equations:
            if eq(values) == 0:
                solved += 1
            else:
                continue
        if solved == len(equations):
            return vector(list(guesses) + constraints)
    return vector([])


def check_solution(equations, solution, reduced=False):
    if reduced:
        for eq in equations:
            assert eq(0, *solution) == 0
    else:
        for eq in equations:
            assert eq(*solution) == 0
    print("The solution is correct")


def count_monomials(equations):
    monomials = set()
    for eq in equations:
        for mon in eq.monomials():
            monomials.add(mon)
    return sorted(list(monomials))


def main():
    q = 4
    m = 4
    n = 8
    uov = UOV(q, m, n)
    verbose = False
    if verbose:
        print("q:", q, ", m:", m, ", n:", n)
        print("Reduced?", uov.reduced)
    equations, matrices = uov.intersection_attack(
        advanced=True, verbose=verbose)
    print("Number of equations:", len(equations))
    print("Number of monomials:", len(count_monomials(equations)))
    if verbose:
        print("")
        print("The system to be solved:")
        for eq in equations:
            print(eq)
        print("")

    solution = guess_solve(equations, q, m, n, xx=uov.xx,
                           advanced=True, reduced=uov.reduced, verbose=verbose)

    if solution == vector([]):
        print("No solution found")
    else:
        print("Solution found:", solution)
        check_solution(equations, solution, uov.reduced)
        for matrix in matrices:
            transformed_solution = matrix * solution
            if uov.reduced:
                transformed_solution = vector([0] + list(transformed_solution))
            if verbose:
                print("Does ", transformed_solution, " lie in O?")
            print(uov.test_oil_space_membership(transformed_solution))


if __name__ == '__main__':
    main()
