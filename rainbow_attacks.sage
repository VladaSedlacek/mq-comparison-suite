from itertools import product
from pathlib import Path
import click
import os


def get_polar_form(Q):
    return Q + Q.transpose()


class Rainbow():
    """A class for the Rainbow scheme."""

    def __init__(self, q, m, n, o2, debug=True):
        assert o2 < m and m < n
        self.debug = debug
        self.q = q
        F = GF(q)
        self.ext_deg = F.modulus().degree()
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
        self.R.inject_variables(verbose=False)
        self.support_ring = PolynomialRing(F, ['y%s' % p for p in range(
            1, n + 1)] + self.support_minors_variables, order="lex")
        self.support_ring.inject_variables(verbose=False)
        self.weil_ring = PolynomialRing(
            F, ['w%s_%s' % (p1, p2) for p1, p2 in product(range(1, n - m - 1), range(self.ext_deg))], order="lex")
        self.weil_ring.inject_variables(verbose=False)
        self.ww = vector(self.weil_ring.gens())
        self.ww_parts = [self.ww[i:i + 4] for i in range(0, 4 * (n - m), 4)]
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
        for _ in range(m - o2):
            Q = Matrix(self.F, n)
            for i in range(n):
                for j in range(n):
                    if i >= n - m or j >= n - o2:
                        continue
                    Q[i, j] = self.F.random_element()
            FF.append(Q)
        for _ in range(m - o2, m):
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
        W_prime_basis = [self.V2([0] * m) for _ in range(o2)]
        for i in range(o2):
            W_prime_basis[i][i + m - o2] = 1
        O1_basis = [self.T.inverse() * o for o in O1_prime_basis]
        O1 = self.V.subspace(O1_basis)
        O2_basis = [self.T.inverse() * o for o in O2_prime_basis]
        O2 = self.V.subspace(O2_basis)
        W_basis = [self.S * o for o in W_prime_basis]
        W = self.V2.subspace(W_basis)

        if self.debug:
            for e in O1_prime_basis:
                assert vector([e * Q * e for Q in self.FF]
                              ) in self.V2.subspace(W_prime_basis)
            for e in O2_prime_basis:
                assert [e * Q * e for Q in self.FF] == [0] * m
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
            assert O1.dimension() == m
            assert O2.dimension() == o2
            assert W.dimension() == o2
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
        max_var_index = n
        if reduce_dimension:
            max_var_index -= o2 - 1
        guessed_vars = [self.F.random_element()
                        for _ in range(n - max_var_index)]
        all_vars = vector(list(yy[:max_var_index]) + guessed_vars)

        def Lx(self, x):
            rows = []
            for e in self.V.basis():
                rows.append([e * M * x for M in self.MM])
            return matrix(self.support_ring, rows)

        Les = [Lx(self, e) for e in self.V.basis()]
        Ly = linear_combination(all_vars, Les)

        if verbose:
            print("Ly:\n", Ly)
        if debug:
            # Check that for random y from O2, the conditions hold (after dimension reduction)
            for _ in range(2 ^ min(len(guessed_vars), 10)):
                y = self.O2.random_element()
                if y[max_var_index:] == all_vars[max_var_index:]:
                    assert Ly(*y, *self.cc).rank() <= o2
                    for row in [row for row in Ly(*y, *self.cc).rows()]:
                        assert vector(self.F(el) for el in row) in self.W

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
            equations.append(all_vars * P * all_vars)

        return equations, guessed_vars


# evaluate Multivariate map and Differential
def Eval(F,x):
    return vector([ x*M*x for M in F])


def Differential(F,x,y):
    return vector([ (x*M*y) + (y*M*x)  for M in F ])


def elt_to_str(a):
    if q == 16:
        return str(hex(sum([2**i * a.polynomial()[i].lift() for i in range(4)])))[2:]
    if q.is_prime() and q < 16:
        return str(hex(a))[2:]
    return str(a)


def UD_to_string(M):
    S = ""
    for i in range(M.ncols()):
        for j in range(i + 1):
            S += elt_to_str(M[j, i]) + ' '
    S += ';\n'
    return S


def differential_attack(self):
    '''Adapted from https://github.com/WardBeullens/BreakingRainbow'''
    q, m, n, o2 = self.q, self.m, self.n, self.o2
    global attempts
    attempts = 0

    # pick a random vector x
    x = vector([F.random_element() for i in range(n)])
    while Eval(self.PP, x)[0] == 0:
        x = vector([F.random_element() for i in range(n)])

    # compute linear map D_x = P'(x,.)
    D_x = Matrix(F, [Differential(self.PP, x, b) for b in (F ^ n).basis])
    D_x_ker = Matrix(D_x.kernel().basis())

    if q % 2 == 0:
        D_x_ker[0] = x

    if D_x_ker.rank() != n - m:
        return Attack(self)

    attempts += 1

    Sol = None
    if not self.O2 is None:
        V = F ^ n
        I = V.span(D_x_ker).intersection(V.span(self.O2.transpose()))
        if I.dimension() == 0:
            print("Attack would fail. resample x")
            return Attack(self)

        print("Intersection has dimension:", I.dimension())
        Sol = I.basis()[0]

        Sol = D_x_ker.transpose().solve_right(Sol)

        if Sol[-1] == 0:
            print("last entry is zero, resample x")
            return Attack(self)

        # scale to solution to have last coordinate zero, reducing dim
        Sol = Sol / Sol[-1]

        print("Good D_x found after %d attempts." % attempts)

        print("The expected solution is:", Sol)
        print("Exp sol in hex format:", [elt_to_str(s) for s in Sol])

    # Compose smaller system D_x(o)= 0 and P(o) = 0
    SS = [D_x_ker * M * D_x_ker.transpose() for M in self.PP]

    for s in SS:
        Make_UD(s)

    if not Sol is None:
        assert Eval(SS, Sol) == vector(m * [0])

    if q % 2 == 0:
        Px = Eval(self.PP, x)
        # define Y by setting first coord to 0?
        SSS = [(SS[i] * Px[0] + SS[0] * Px[i])[1:, 1:]
               for i in range(1, len(SS))]

        if Sol is not None:
            assert Eval(SSS, Sol[1:]) == vector((m - 1) * [0])

        SS_orig = SS
        SS = SSS

        # In the real attack, YSol is found externally via a solver
        YSol = vector([0] + list(Sol[1:]))
        alpha = Eval([SS_orig[0] / Px[0]], YSol)[0]

        xt = D_x_ker.transpose().solve_right(x)
        assert xt == vector([1] + (n - m - 1) * [0])
        assert Eval(SS_orig, xt) == Px
        assert Eval(SS_orig, YSol) == alpha * Px

        NewSol = xt + 1 / alpha.sqrt() * YSol
        NewSol = NewSol / NewSol[-1]
        assert NewSol == Sol
        assert Eval(SS_orig, NewSol) == vector(m * [0])
        print("New sol in hex format:", [elt_to_str(s) for s in NewSol])

    fname = 'system_' + str(len(SS)) + '-' + str(SS[0].ncols() - 1) + '.xl'
    f = open(fname, 'w')
    for s in SS:
        f.write(UD_to_string(s))
    f.close()
    print("System written to: " + fname)
    print("Use block Wiedemann XL algorithm of Niederhagen to find a solution:")
    print("http://polycephaly.org/projects/xl")
    xl_path = Path("..", "xl")
    make_command = "make -C {} Q={} M={} N={}".format(
        str(xl_path), str(q), str(len(SS)), str(SS[0].ncols() - 1))
    os.system(make_command)

    z = F.gens()[0]
    zs = [z ^ i for i in range(e)]
    yy_weil = vector([linear_combination(w, zs) for w in self.ww_parts])
    yy_weil_affine = vector(list(yy_weil[1:-1]) + [1])
    equations = [yy_weil_affine * s * yy_weil_affine for s in SS]
    assert len(equations) == m - 1
    equations_weil = [weil_decomposition(eq) for eq in equations]
    equations_final = [delete_powers(w_eq)
                       for eqs in equations_weil for w_eq in eqs]

    save_system(equations_final,
                "system_{}-{}.mq".format(m - 1, n - m - 2))


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
    oil_vector = vector(list(solution[: max_nonzero_index]
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


def delete_powers(eq):
    return sum([radical(mon) for mon in eq.monomials()])


def save_system(rainbow, equations, guessed_vars, reduce_dimension, file_path):
    var_set = set()
    for eq in equations:
        if eq == 0:
            continue
        for var in eq.variables():
            var_set.add(var)
    var_list = [str(var) for var in sorted(var_set)[:: -1]]
    variables = ', '.join(var_list)
    max_var_index = rainbow.n
    if reduce_dimension:
        max_var_index -= rainbow.o2 - 1
    guessed = ', '.join(["{0}={1}".format(var, value) for var, value in zip(
        rainbow.yy[max_var_index:], guessed_vars)])

    with open(file_path, 'w') as file:
        file.write("# Variables:\n")
        file.write(variables + "\n\n")
        file.write("# Guessed variables:\n")
        file.write("# " + guessed + "\n\n")
        file.write("# Equations:\n")
        for eq in equations:
            file.write(str(eq) + "\n")


def save_setup(rainbow, setup_path):
    with open(setup_path, 'w') as file:
        file.write("# O1:" + str(rainbow.O1) + "\n\n")
        file.write("# O2:" + str(rainbow.O2) + "\n\n")
        file.write("# W:" + str(rainbow.W) + "\n\n")


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


def try_toy_solution(rainbow, equations, attack_type, reduce_dimension):
    print("Constructing toy solution...")
    if attack_type == 'minrank':
        solution = vector(list(rainbow.O2.random_element()) +
                          [0] * len(rainbow.support_minors_indices))
        success = check_rectangular_minrank_attack_success(
            equations, solution, rainbow, reduce_dimension=reduce_dimension)
    elif attack_type == 'intersection':
        solution = vector([0] * n + list(rainbow.W.complement().basis()[0]))
        success = check_intersection_attack_success(
            equations, solution, rainbow)
    if success:
        print("Attack successful!")
    else:
        print("Attack not successful :(")


@ click.command()
@ click.option('--q', default=2, help='the field order', type=int)
@ click.option('--o2', default=2, help='the oil subspace dimension', type=int)
@ click.option('--m', default=4, help='the number of equations', type=int)
@ click.option('--n', default=8, help='the number of variables', type=int)
@ click.option('-n', '--no_solve', default=False, is_flag=True, help='try to solve the system or not')
@ click.option('--mq_path', default=Path("..", "mq"), help='the path the MQ solver: https://gitlab.lip6.fr/almasty/mq', type=str)
@ click.option('-h', '--inner_hybridation', default="-1", help='the number of variable that are not guessed', type=int)
@ click.option('-v', '--verbose', default=False, is_flag=True, help='control the output verbosity')
@ click.option('-r', '--reduce_dimension', default=True, is_flag=True, help='reduce the dimension when possible')
@ click.option('-w', '--weil_descent', default=True, is_flag=True, help='use Weil descent when possible')
@ click.option('-t', '--attack_type', default='minrank', type=click.Choice(['minrank', 'intersection'], case_sensitive=False), help='use either the rectangular MinRank attack or the intersection attack')
def main(q, o2, m, n, no_solve, mq_path, inner_hybridation, verbose, reduce_dimension, weil_descent, attack_type):
    boolean = q % 2 == 0
    set_random_seed(0)
    rainbow = Rainbow(q, m, n, o2)
    if verbose:
        print("O1:", rainbow.O1, "\n")
        print("O2:", rainbow.O2, "\n")
        print("W:", rainbow.W, "\n")

    if attack_type == 'minrank':
        if verbose:
            print("Mounting the rectangular MinRank attack...")
        equations, guessed_vars = rainbow.rectangular_minrank_attack(
            reduce_dimension=reduce_dimension, verbose=verbose)
    elif attack_type == 'intersection':
        if verbose:
            print("Mounting the intersection_attack...")
        equations, _, matrices = rainbow.intersection_attack()
        guessed_vars = []

    if weil_descent:
        if verbose:
            print("\nPerforming Weil descent...")
        equations = [eq_desc for equations_desc in [weil_decomposition(
            eq) for eq in equations] for eq_desc in equations_desc]
    if boolean:
        assert weil_descent or is_prime(q)
        equations = [delete_powers(eq) for eq in equations]
    if verbose:
        print("Number of equations:", len(equations))
        print("Number of monomials:", len(count_monomials(equations)))
        print("")
        print("The system to be solved:")
        for eq in equations:
            print(eq)
        print("")

    eq_path = Path(
        'systems', "rainbow_q_{0}_o2_{1}_m_{2}_n_{3}.eq".format(q, o2, m, n))
    setup_path = Path(
        'systems', "rainbow_q_{0}_o2_{1}_m_{2}_n_{3}.stp".format(q, o2, m, n))
    save_system(rainbow, equations, guessed_vars, reduce_dimension, eq_path)
    save_setup(rainbow, setup_path)
    print("Saving the equation system into", eq_path)

    if not no_solve:
        print("Starting the MQ solver")
        if inner_hybridation == -1:
            inner_hybridation_arg = ""
        else:
            inner_hybridation_arg = " --inner-hybridation " + \
                str(inner_hybridation)
        mq_solve_command = "{}{} < {}".format(
            str(Path(mq_path, "monica_vector")), inner_hybridation_arg, str(eq_path))
        os.system(mq_solve_command)


if __name__ == '__main__':
    main()
