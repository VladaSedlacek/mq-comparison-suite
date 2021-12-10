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
        self.k = find_max_k(self.m, self.n, verbose=False)
        self.reduced = self.q % 2 == 0 and self.n % 2 == 1
        self.V = VectorSpace(F, n)
        self.V2 = VectorSpace(F, m)
        self.R = PolynomialRing(F, ['x%s' % p for p in range(
            1, n + 1)], order="neglex")
        self.R.inject_variables()
        self.xx = vector(self.R.gens()[:n])
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
        O_prime_basis = [self.V([0] * n) for _ in range(m)]
        for i in range(m):
            O_prime_basis[i][i + n - m] = 1
        O_basis = [self.T.inverse() * o for o in O_prime_basis]
        O = self.V.subspace(O_basis)

        if self.debug:
            for o in O_prime_basis:
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

    def intersection_attack(self, verbose=False):
        m, n = self.m, self.n
        assert 2 * m <= n and n < 3 * m
        if self.reduced:
            m -= 1
            n -= 1
            MM = [M[1:, 1:] for M in self.MM]
            PP = [P[1:, 1:] for P in self.PP]
            xx = (self.xx)[1:]
        else:
            MM = self.MM
            PP = self.PP
            xx = self.xx
        LL = []
        combinations = []
        for i in range(self.k):
            while True:
                coefficients = self.V2.random_element()
                L = linear_combination(coefficients, MM)
                if L.is_invertible():
                    LL.append(L)
                    combinations.append(coefficients)
                    break
        equations = []
        redundant = []
        for i in range(self.k):
            u = LL[i].inverse() * xx
            for P in PP:
                equations.append(u * P * u)
            for j in range(i + 1, self.k):
                v = LL[j].inverse() * xx
                for l, M in enumerate(MM):
                    eq = u * M * v
                    nonzero_index_i, nonzero_index_j = first_different_nonzero_indices(
                        combinations[i], combinations[j])
                    if l != nonzero_index_i and l != nonzero_index_j:
                        equations.append(eq)
                    else:
                        redundant.append(eq)
        matrices = [L.inverse() for L in LL]
        return equations, redundant, matrices


def first_nonzero_index(it):
    for i, _ in enumerate(it):
        if it[i] != 0:
            return i
    return None


def first_different_nonzero_indices(it1, it2):
    i = first_nonzero_index(it1)
    j = first_nonzero_index(it2)
    if i == j:
        j = first_nonzero_index(
            [0] + list(it2))
    assert i != j
    return i, j


def linear_combination(coefficients, objects):
    assert len(coefficients) == len(objects)
    return sum(c * o for c, o in zip(coefficients, objects))


def find_max_k(m, n, verbose=False):
    if n == 2 * m:
        return ceil(sqrt(m))
    k = 2
    while True:
        if verbose:
            print("current k:", k, ", n/m:", (n / m).numerical_approx(digits=3),
                  ", (2 * k - 1) / (k - 1):", ((2 * k - 1) / (k - 1)).numerical_approx(digits=3))
        if n >= (2 * k - 1) / (k - 1) * m:
            k -= 1
            break
        k += 1
    if verbose:
        print("k:", k)
    assert n > k * m - (k - 1) * (n - m)
    return k


def guess_solve(equations, uov, verbose=False):
    q, m, n, k, xx, reduced = uov.q, uov.m, uov.n, uov.k, uov.xx, uov.reduced
    total_count = ZZ(m * k * (k + 1) / 2)
    redundant_count = k * (k - 1)
    assert len(equations) == total_count - redundant_count
    tail = k * m - (k - 1) * (n - m)
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
    n = 9
    uov = UOV(q, m, n)
    k = uov.k
    verbose = False
    if verbose:
        print("q:", q, ", m:", m, ", n:", n, " k:", k)
        print("Reduced?", uov.reduced)
    equations, _, matrices = uov.intersection_attack(verbose=verbose)
    print("Number of equations:", len(equations))
    print("Number of monomials:", len(count_monomials(equations)))
    if verbose:
        print("")
        print("The system to be solved:")
        for eq in equations:
            print(eq)
        print("")

    solution = guess_solve(equations, uov, verbose=verbose)

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
