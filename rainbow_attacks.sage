from itertools import product


def get_polar_form(Q):
    return Q + Q.transpose()


class Rainbow():
    """A class for the Rainbow scheme."""

    def __init__(self, q, m, n, o2, debug=True):
        assert o2 < m and m < n
        self.debug = debug
        self.q = q
        F = GF(q)
        self.F = F
        self.m = m
        self.n = n
        self.o2 = o2
        self.k = find_max_k(m=self.o2, n=self.n, verbose=False)
        self.reduced = self.q % 2 == 0 and self.n % 2 == 1
        self.V = VectorSpace(F, n)
        self.V2 = VectorSpace(F, m)
        self.support_minors_indices = [
            str(c) for c in Combinations([1..self.m], self.o2)]
        self.support_minors_variables = [
            'c%s' % s for s in [1..len(list(self.support_minors_indices))]]
        self.support_minors_dict = dict(
            zip(self.support_minors_indices, self.support_minors_variables))
        self.R = PolynomialRing(F, ['x%s' % p for p in range(
            1, n + 1)] + ['v%s' % p for p in range(1, m + 1)], order="lex")
        self.R.inject_variables()
        self.support_ring = PolynomialRing(F, ['y%s' % p for p in range(
            1, n + 1)] + self.support_minors_variables, order="lex")
        self.support_ring.inject_variables()
        self.xx = vector(self.R.gens()[: n])
        self.vv = vector(self.R.gens()[n:])
        self.yy = vector(self.support_ring.gens()[: n])
        self.cc = vector(self.support_ring.gens()[n:])
        self.FF = self.construct_central_map()
        self.T, self.S, self.PP, self.MM = self.hide_central_map()
        self.O1, self.O2, self.W = self.find_subspaces()

    def construct_central_map(self):
        m, n, o2 = self.m, self.n, self.o2
        FF = []
        for _ in range(o2):
            Q = Matrix(self.F, n)
            for i in range(n):
                for j in range(n):
                    if i >= n - o2 or j >= n - m:
                        continue
                    Q[i, j] = self.F.random_element()
            FF.append(Q)
        for _ in range(o2, m):
            Q = Matrix(self.F, n)
            for i in range(n):
                for j in range(n):
                    if j >= n - o2:
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
        while True:
            S = random_matrix(self.F, m)
            if S.rank() == m:
                break
        PP_untransformed = [T.transpose() * Q * T for Q in self.FF]
        PP = PP_untransformed
        PP = [linear_combination(S_row, PP_untransformed)
              for S_row in S.rows()]
        MM = [get_polar_form(P) for P in PP]
        return T, S, PP, MM

    def find_subspaces(self):
        m, n, o2 = self.m, self.n, self.o2
        O1_prime_basis = [self.V([0] * n) for _ in range(m)]
        for i in range(m):
            O1_prime_basis[i][i + n - m] = 1
        O2_prime_basis = [self.V([0] * n) for _ in range(o2)]
        for i in range(o2):
            O2_prime_basis[i][i + n - o2] = 1
        W_prime_basis = [self.V2([0] * m) for _ in range(m - o2)]
        for i in range(m - o2):
            W_prime_basis[i][i + o2] = 1
        O1_basis = [self.T.inverse() * o for o in O1_prime_basis]
        O1 = self.V.subspace(O1_basis)
        O2_basis = [self.T.inverse() * o for o in O2_prime_basis]
        O2 = self.V.subspace(O2_basis)
        W_basis = [self.S * o for o in W_prime_basis]
        W = self.V2.subspace(W_basis)

        if self.debug:
            for e in O1_basis:
                assert vector([e * P * e for P in self.PP]
                              ) in self.V2.subspace(W_basis)
            for e in O2_basis:
                assert [e * P * e for P in self.PP] == [0] * m
                for f in self.V.basis():
                    assert vector([f * M * e for M in self.MM]
                                  ) in self.V2.subspace(W_basis)
            for v in W.complement().basis():
                Mv = linear_combination(v, self.MM)
                assert O2.is_subspace(Mv.kernel())
        return O1, O2, W

    def intersection_attack(self, verbose=False):
        if self.reduced:
            print("Reduced variant not implemented yet")
            exit()
        MM = self.MM
        PP = self.PP
        xx = self.xx
        vv = self.vv
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
            for e in self.V.basis():
                eq = linear_combination(vv, [u * M * e for M in MM])
                equations.append(eq)
        matrices = [L.inverse() for L in LL]
        return equations, redundant, matrices

    def rectangular_minrank_attack(self, reduce_dimension=False, debug=True, verbose=True):
        m, n, o2 = self.m, self.n, self.o2
        yy = self.yy
        max_nonzero_index = n
        if reduce_dimension:
            max_nonzero_index -= o2 - 1

        def Lx(self, x):
            rows = []
            for e in self.V.basis():
                rows.append([e * M * x for M in self.MM])
            return matrix(rows)

        Les = [Lx(self, e) for e in self.V.basis()[:max_nonzero_index]]

        Ly = matrix(self.support_ring, n, m)
        for i in range(max_nonzero_index):
            summand_matrix_rows = []
            for j, row in enumerate(Les[i].rows()):
                summand_matrix_rows.append([yy[i] * el for el in row])
            Ly += matrix(summand_matrix_rows)

        if debug:
            for _ in range(10):
                y = self.O2.random_element()
                assert Ly(*y, *self.cc).rank() <= o2

        if verbose:
            print("Ly:\n", Ly)

        equations = []

        # Add equations for Support Minors Modeling.
        for j in range(n):
            rj = Ly[j]
            for columns in Combinations([1..self.m], self.o2 + 1):
                eq = 0
                for main_index in columns:
                    support_indices = [c for c in columns if c != main_index]
                    cs = self.support_ring(
                        self.support_minors_dict[str(support_indices)])
                    eq += (-1) ^ (main_index - 1) * rj[main_index - 1] * cs
                equations.append(eq)

        # Add quadratic equations for oil subspace membership.
        for P in self.PP:
            equations.append(yy[:max_nonzero_index] * P[:max_nonzero_index,
                                                        :max_nonzero_index] * yy[:max_nonzero_index])

        return equations


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
    if n == 3 * m:
        return 2
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


def check_solution(equations, solution, reduced=False):
    if len(solution) == 0:
        print("No solution found")
        return False
    if reduced:
        solution = vector([0] + list(solution))
    for eq in equations:
        assert eq(*solution) == 0
    print("The solution is correct")
    return True


def check_intersection_attack_success(equations, solution, rainbow, verbose=True):
    print("Solution found:", solution)
    success = check_solution(equations, solution, rainbow.reduced)
    v = solution[rainbow.n:]
    Mv = linear_combination(v, rainbow.MM)
    return success and rainbow.O2.is_subspace(Mv.kernel())


def check_rectangular_minrank_attack_success(equations, solution, rainbow, reduce_dimension=False, verbose=True):
    print("Solution found:", solution)
    success = check_solution(equations, solution, rainbow.reduced)
    n = rainbow.n
    o2 = rainbow.o2
    max_nonzero_index = n
    if reduce_dimension:
        max_nonzero_index -= o2 - 1
    oil_vector = vector(list(solution[:max_nonzero_index]
                             ) + [0] * (n - max_nonzero_index))
    return success and oil_vector in rainbow.O2


def count_monomials(equations):
    monomials = set()
    for eq in equations:
        for mon in eq.monomials():
            monomials.add(mon)
    return sorted(list(monomials))


def compare_variants(n, m, o2):
    print("m = {0}, n = {1}, o2 = {2}".format(m, n, o2))
    smm_equations = n * binomial(m, o2 + 1) + m
    smm_variables = binomial(m, o2) + n - o2 + 1
    alt_equations = m * (n + 1)
    alt_variables = m * (n + o2) + n - o2 + 1
    print("\tSMM: {0} equations, {1} variables \n\tratio: {2}, difference: {3}".format(
        smm_equations, smm_variables, (smm_variables / smm_equations).numerical_approx(digits=3), smm_variables - smm_equations))  # %
    print("\tALT: {0} equations, {1} variables \n\tratio: {2}, difference: {3}".format(
        alt_equations, alt_variables, (alt_variables / alt_equations).numerical_approx(digits=3), alt_variables - alt_equations))


def print_nist_comparisons():
    # Rainbow notation mapping:
    # v1 + o2 + o1 -> n
    #      o2 + o1 -> m
    #           o1 -> o2
    print("\n\nNIST parameters over F_16 (Ia, IVa, VIa):")
    compare_variants(n=96, m=64, o2=32)
    compare_variants(n=152, m=96, o2=48)
    compare_variants(n=204, m=128, o2=64)
    print("\n\nNIST parameters over F_31 (Ib, IIIb, VIb:")
    compare_variants(n=92, m=56, o2=28)
    compare_variants(n=144, m=80, o2=32)
    compare_variants(n=196, m=112, o2=56)
    print("\n\nNIST parameters over F_256 (Ic, IIIc, Vc):")
    compare_variants(n=88, m=48, o2=24)
    compare_variants(n=140, m=72, o2=36)
    compare_variants(n=188, m=96, o2=48)
    print("\n\nNew suggested parameters over F_16 (I):")
    compare_variants(m=68, n=109, o2=36)


def main():
    q = 2
    o2 = 2
    m = 4
    n = 8
    verbose = False
    minrank = True
    reduce_dimension = False
    rainbow = Rainbow(q, m, n, o2)
    if verbose:
        print("O1:", rainbow.O1, "\n")
        print("O2:", rainbow.O2, "\n")
        print("W:", rainbow.W, "\n")

    if minrank:
        print("Mounting the rectangular MinRank attack...")
        equations = rainbow.rectangular_minrank_attack(
            reduce_dimension=reduce_dimension, verbose=False)
    else:
        print("Mounting the intersection_attack...")
        equations, _, matrices = rainbow.intersection_attack()
    print("Number of equations:", len(equations))
    print("Number of monomials:", len(count_monomials(equations)))
    print("")
    print("The system to be solved:")
    for eq in equations:
        print(eq)
    print("")

    print("Constructing toy solution...")
    if minrank:
        solution = vector(list(rainbow.O2.random_element()) +
                          [0] * len(rainbow.support_minors_indices))
        success = check_rectangular_minrank_attack_success(
            equations, solution, rainbow, reduce_dimension=reduce_dimension)
    else:
        solution = vector([0] * n + list(rainbow.W.complement().basis()[0]))
        success = check_intersection_attack_success(
            equations, solution, rainbow)
    if success:
        print("Attack successful!")
    else:
        print("Attack not successful :(")

    print_nist_comparisons()


if __name__ == '__main__':
    main()
