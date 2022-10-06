#!/usr/bin/env sage

from itertools import product
from pathlib import Path
import click
import re
from invoke_solver import invoke_solver
from utils import get_eq_format


def get_polar_form(Q):
    return Q + Q.transpose()


class Rainbow():
    """A class for the Rainbow scheme."""

    def __init__(self, q, m, n, o2, seed=0, debug=True, support=False):
        assert o2 < m and m < n
        self.seed = seed
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
        self.order = "lex"
        self.R = PolynomialRing(F, ['x%s' % p for p in range(
            1, n + 1)] + ['v%s' % p for p in range(1, m + 1)], order=self.order)
        self.R.inject_variables(verbose=False)
        self.weil_ring = PolynomialRing(
            F, ['w%s_%s' % (p1, p2) for p1, p2 in product(range(1, n - m + 1), range(self.ext_deg))], order='degrevlex')
        self.weil_ring.inject_variables(verbose=False)
        self.ww = vector(self.weil_ring.gens())
        self.ww_parts = [self.ww[i:i + self.ext_deg]
                         for i in range(0, self.ext_deg * (n - m), self.ext_deg)]
        self.xx = vector(self.R.gens()[: n])
        self.vv = vector(self.R.gens()[n:])
        if support:
            self.support_minors_indices = [
                str(c) for c in Combinations([1..self.m], self.o2)]
            self.support_minors_variables = [
                'c%s' % s for s in [1..len(list(self.support_minors_indices))]]
            self.support_minors_dict = dict(
                zip(self.support_minors_indices, self.support_minors_variables))
            self.support_ring = PolynomialRing(F, ['y%s' % p for p in range(
                1, n + 1)] + self.support_minors_variables, order=self.order)
            self.support_ring.inject_variables(verbose=False)
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

    def differential_attack(self, debug=False, verbose=False):
        '''Adapted from https://github.com/WardBeullens/BreakingRainbow'''
        q, m, n, o2 = self.q, self.m, self.n, self.o2
        global attempts
        attempts = 0

        # pick a random vector x
        x = vector([self.F.random_element() for i in range(n)])
        while Eval(self.PP, x)[0] == 0:
            x = vector([self.F.random_element() for i in range(n)])

        # compute linear map D_x = P'(x,.)
        D_x = Matrix(self.F, [Differential(self.PP, x, b)
                              for b in (self.F ^ n).basis()])
        D_x_ker = Matrix(D_x.kernel().basis())

        if q % 2 == 0:
            D_x_ker[0] = x

        if D_x_ker.rank() != n - m:
            return self.differential_attack()

        attempts += 1

        Sol = None
        if not self.O2 is None:
            V = self.F ^ n
            I = V.span(D_x_ker).intersection(self.O2)
            if I.dimension() == 0:
                if verbose:
                    print("\tAttack would fail. Resampling x...")
                return self.differential_attack()
            if verbose:
                print("\tIntersection has dimension:", I.dimension())
            Sol = I.basis()[0]

            Sol = D_x_ker.transpose().solve_right(Sol)

            if Sol[-1] == 0:
                if verbose:
                    print("\tLast entry is zero, resampling x...")
                return self.differential_attack()

            # scale to solution to have last coordinate zero, reducing dim
            Sol = Sol / Sol[-1]
            if verbose:
                print("\tGood D_x found after %d attempts." % attempts)

                print("\tThe expected solution is:", Sol)
                print("\tIn hex format:", [elt_to_str(q, s) for s in Sol])

        # Compose smaller system D_x(o)= 0 and P(o) = 0
        SS = [D_x_ker * M * D_x_ker.transpose() for M in self.PP]

        for s in SS:
            Make_UD(s)

        if not Sol is None:
            assert Eval(SS, Sol) == vector(m * [0])

        weil_coeff_list = []
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

            if alpha == 0:
                NewSol = xt
            else:
                NewSol = xt + 1 / alpha.sqrt() * YSol
            if NewSol[-1] != 0:
                NewSol = NewSol / NewSol[-1]
            # assert NewSol == Sol
            # assert Eval(SS_orig, NewSol) == vector(m * [0])

            # Perform the Weil descent
            z = self.F.gens()[0]
            zs = [z ^ i for i in range(self.ext_deg)]
            yy_weil = vector([linear_combination(w, zs)
                              for w in self.ww_parts])
            yy_weil_affine = vector(list(yy_weil[1:-1]) + [1])
            equations = [yy_weil_affine * s * yy_weil_affine for s in SS]
            assert len(equations) == m - 1
            equations_weil = [
                w_eq for eq in equations for w_eq in weil_decomposition(eq)]
            # intentionally do not call delete_powers at this point to be consistent with XL representation
            equations_final = equations_weil

            # Get Weil coefficients from the equations
            weil_vars = [
                var for part in self.ww_parts for var in part][self.ext_deg:-self.ext_deg]
            # get all monomials that are at most quadratic, in grevdeglex order
            weil_vars.append(1)
            l = len(weil_vars)
            rows = [i * [0] + (l - i) * [1] for i in range(l)]
            upper_matrix = Matrix(self.F, rows)
            max_quadratic = vector(weil_vars) * \
                upper_matrix * vector(weil_vars)
            for eq in equations_final:
                weil_coeffs = []
                for mon in max_quadratic.monomials():
                    weil_coeffs.append(eq.monomial_coefficient(mon))
                assert linear_combination(
                    max_quadratic.monomials(), weil_coeffs) == eq
                weil_coeff_list.append(weil_coeffs)
        return SS, equations_final, weil_coeff_list, Sol


# evaluate Multivariate map and Differential
def Eval(F, x):
    return vector([x * M * x for M in F])


def Differential(F, x, y):
    return vector([(x * M * y) + (y * M * x) for M in F])


# makes a matrix upper diagonal
def Make_UD(M):
    n = M.ncols()
    K = M[0, 0].parent()
    for i in range(n):
        for j in range(i + 1, n):
            M[i, j] += M[j, i]
            M[j, i] = K(0)


def elt_to_str(q, a):
    if q == 16:
        return str(hex(sum([2**i * a.polynomial()[i].lift() for i in range(4)])))[2:]
    if ZZ(q).is_prime() and q < 16:
        return str(hex(a))[2:]
    return str(a)


def str_to_elt(q, s):
    return GF(q).fetch_int(int(s, q))


def four_bits_to_str(bit_list):
    return str(hex(int(''.join(map(str, bit_list)).encode(), 2)))[2:]


def UD_to_string(q, M):
    # this corresponds to the degrevlex order
    S = ""
    for i in range(M.ncols()):
        for j in range(i + 1):
            S += elt_to_str(q, M[j, i]) + ' '
    S += ';\n'
    return S


def weil_coeff_list_to_string(weil_coeff_list, deg):
    S = ""
    for weil_coeffs in weil_coeff_list:
        weil_coeffs += ((-len(weil_coeffs)) % deg) * [0]
        assert len(weil_coeffs) % deg == 0
        for i in range(len(weil_coeffs) / deg):
            bit_list = weil_coeffs[deg * i:deg * i + deg]
            S += four_bits_to_str(bit_list) + ' '
        S += ';\n'
    return S


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


def equations_to_degrevlex_str(equations):
    q = equations[0].base_ring().order()
    var_set = set().union(*[eq.variables() for eq in equations])
    var_list = sorted(var_set)[:: -1] + [1]
    degrevlex_mons = []
    # degrevlex corresponds to reading the upper part of the quadratic form matrix column-wise
    for j in range(len(var_list)):
        for i in range(j + 1):
            degrevlex_mons.append(var_list[i] * var_list[j])

    coeff_str = ""
    for eq in equations:
        coeffs = [eq.monomial_coefficient(mon) for mon in degrevlex_mons[:-1]] + [eq.constant_coefficient()]
        coeff_str += " ".join([elt_to_str(q, coeff) for coeff in coeffs]) + " ;\n"
    return coeff_str


def save_system(file_format, file_path, rainbow, equations, verbose=False):
    if file_path.is_file() and verbose:
        print("The file {} already exists!".format(str(file_path)))
        return

    equations = [eq for eq in equations if eq != 0]

    if file_format == 'cb_gpu':
        '''The format for the GPU F2 crossbred solver of Niederhagen, Ning and Yang: https://github.com/kcning/mqsolver/'''
        var_set = set().union(*[eq.variables() for eq in equations if eq != 0])
        M = len(equations)
        N = len(var_set)
        with open(file_path, 'w') as f:
            f.write(
                """Galois Field : GF(2)
Number of variables (n) : {}
Number of polynomials (m) : {}
Seed : {}
Order : graded reverse lex order

*********************\n""".format(N, M, rainbow.seed))
            f.write(equations_to_degrevlex_str(equations))

    elif file_format == 'cb_orig':
        var_set = set().union(*[eq.variables() for eq in equations if eq != 0])
        M = len(equations)
        N = len(var_set)
        with open(file_path, 'w') as f:
            for eq in equations:
                eq_repr = []
                for var_tuple, coeff in eq.dict().items():
                    if coeff == 1:
                        # create an integer whose binary representation corresponds to variables present in the monomial; divide by 2 to disregard the last variable, which should not be present
                        mon_repr = ZZ(''.join([str(k) for k in var_tuple]), 2)
                        assert mon_repr % 2 == 0
                        eq_repr.append(mon_repr / 2)
                for mon_repr in sorted(eq_repr, reverse=True):
                    f.write(str(mon_repr) + "\n")
                f.write(str(-1) + "\n")

    elif file_format == 'cnf':
        var_set = set().union(*[eq.variables() for eq in equations])
        var_list = sorted(var_set)[:: -1]
        var_prod_list = []
        M = len(equations)
        N = len(var_set)
        with open(file_path, 'w') as f:
            f.write("p cnf {} {}\n".format(
                N + binomial(N, 2) + 1, M + 3 * binomial(N, 2) + 1))
            # introduce the constant variable
            f.write("{} 0\n".format(N + binomial(N, 2) + 1))
            # convert ANDs to ORs by introducing new variables
            prod_index = 0
            for i, _ in enumerate(var_list):
                for j in range(i + 1, len(var_list)):
                    prod_index += 1
                    var_prod_list.append(var_list[i] * var_list[j])
                    f.write("{} -{} 0\n".format(i + 1, N + prod_index))
                    f.write("{} -{} 0\n".format(j + 1, N + prod_index))
                    f.write("-{} -{} {} 0\n".format(i +
                                                    1, j + 1, N + prod_index))
            for eq in equations:
                const_present = False
                cnf_line = "x "
                for mon in eq.monomials():
                    # add linear monomials
                    if mon in var_list:
                        cnf_line += "{} ".format(var_list.index(mon) + 1)
                    # add quadratic monomials
                    elif mon in var_prod_list:
                        cnf_line += "{} ".format(
                            N + var_prod_list.index(mon) + 1)
                    else:
                        assert mon == 1
                        const_present = True
                if not const_present:
                    # the right hand side of the equation must correspond to True
                    cnf_line += "{} ".format(str(N + binomial(N, 2) + 1))
                cnf_line += "0\n"
                f.write(cnf_line)

    elif file_format == 'magma':
        q = rainbow.q
        if q == 2:
            # for optimized performance
            ring_type = "BooleanPolynomialRing("
        else:
            ring_type = "PolynomialRing(F, "
        var_set = set().union(*[eq.variables() for eq in equations])
        var_list = [str(var) for var in sorted(var_set)[:: -1]]
        with open(file_path, 'w') as f:
            f.write(f"F := GaloisField({q});\n")
            f.write(f"R<{', '.join(var_list)}> := {ring_type}{len(var_list)}, \"grevlex\");\n")
            f.write(f"I := ideal<R |\n")
            # add field equations
            for v in var_list:
                f.write(f"{v}^{q} - {v},\n")
            # add system equations
            for i, eq in enumerate(equations):
                f.write(f"{eq}")
                if i != len(equations) - 1:
                    f.write(",\n")
                else:
                    f.write("\n")
            f.write(f">;\n")
            # use the F4 algorithm
            f.write("GroebnerBasis(I: Faugere:=true);\n")
            f.write("Variety(I);")

    elif file_format == 'mq':
        var_set = set().union(*[eq.variables() for eq in equations])
        var_list = [str(var) for var in sorted(var_set)[:: -1]]
        variables = ', '.join(var_list)
        max_var_index = rainbow.n
        if reduce_dimension:
            max_var_index -= rainbow.o2 - 1
        if guessed_vars == []:
            guessed = ""
        else:
            guessed = ', '.join(["{0}={1}".format(var, value) for var, value in zip(
                rainbow.yy[max_var_index:], guessed_vars)])

        with open(file_path, 'w') as f:
            f.write("# Variables:\n")
            f.write(variables + "\n#\n")
            f.write("# Guessed variables:\n")
            f.write("# " + guessed + "\n#\n")
            f.write("# Equations:\n")
            for eq in equations:
                f.write(str(eq) + "\n")
        if verbose:
            print("Number of equations:", len(equations))
            print("Number of monomials:", len(
                count_monomials(equations)))

    elif file_format == 'mq_compact':
        with open(file_path, 'w') as f:
            f.write(weil_coeff_list_to_string(
                weil_coeff_list, rainbow.ext_deg))

    elif file_format == 'wdsat':
        var_set = set().union(*[eq.variables() for eq in equations if eq != 0])
        var_list = sorted(var_set)[:: -1]
        var_prod_dict = {v1 * v2: sorted([i + 1, j + 1]) for i, v1 in enumerate(
            var_list) for j, v2 in enumerate(var_list) if v1 != v2}
        M = len(equations)
        N = len(var_set)
        with open(file_path, 'w') as f:
            f.write("p anf {} {}\n".format(N, M))
            for eq in equations:
                const_present = False
                anf_line = "x "
                for mon in eq.monomials():
                    if mon in var_list:
                        anf_line += "{} ".format(var_list.index(mon) + 1)
                    elif mon in var_prod_dict.keys():
                        anf_line += ".2 {} {} ".format(
                            var_prod_dict[mon][0], var_prod_dict[mon][1])
                    else:
                        assert mon == 1
                        const_present = True
                if not const_present:
                    # the right hand side of the equation must correspond to True
                    anf_line += "T "
                anf_line += "0\n"
                f.write(anf_line)

    elif file_format == 'xl':
        '''The format for the block Wiedemann XL solver of Niederhagen: http://polycephaly.org/projects/xl'''
        with open(file_path, 'w') as f:
            f.write(equations_to_degrevlex_str(equations))

    assert file_format in ['cb_gpu', 'cb_orig', 'cnf', 'magma', 'mq', 'mq_compact', 'wdsat', 'xl']
    if verbose:
        print("Equation system written to: " + str(file_path))


def save_setup(rainbow, setup_path, verbose=False):
    if setup_path.is_file() and verbose:
        print("The file {} already exists!".format(str(setup_path)))
        return
    with open(setup_path, 'w') as f:
        f.write("# seed:" + str(rainbow.seed) + "\n")
        f.write("# O1:" + str(rainbow.O1) + "\n\n")
        f.write("# O2:" + str(rainbow.O2) + "\n\n")
        f.write("# W:" + str(rainbow.W) + "\n\n")


def save_solution(solution, solution_path, verbose=False):
    if solution_path.is_file() and verbose:
        print("The file {} already exists!".format(str(solution_path)))
        return
    with open(solution_path, 'w') as f:
        for s in solution:
            f.write(str(s) + "\n")


def load_solution(solution_path, q):
    try:
        F = GF(q)
        with open(solution_path, 'r') as f:
            lines = f.readlines()
            solution = [F(s) for s in lines]
            return vector(solution)
    except Exception as e:
        print("An error ocurred during loading the solution: ", e)
        exit()


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


def compute_system_size(q, m, n, o2, attack_type):
    '''Return the number of equations and the number of variables'''
    if attack_type == 'differential':
        if q % 2 == 0:
            return m - 1, n - m - 2
        return m, n - m - 1
    if attack_type == 'minrank':
        # TODO: this should be different for even q
        return n * binomial(m, o2 + 1) + m, n - o2 + 1 + binomial(m, o2)
    if attack_type == 'intersection':
        k = find_max_k(o2, n)
        return binomial(k + 1, 2) * m - k * (k - 1), n + m
    return None, None


def mount_attack(rainbow, attack_type, M, N, reduce_dimension=False, verbose=False):
    SS = []
    guessed_vars = []
    weil_coeff_list = []
    if attack_type == 'differential':
        if verbose:
            print("Mounting the differential attack...")
        SS, equations, weil_coeff_list, solution = rainbow.differential_attack(
            debug=False, verbose=verbose)
        assert M == len(SS)
        assert N == SS[0].ncols() - 1
    elif attack_type == 'minrank':
        if verbose:
            print("Mounting the rectangular MinRank attack...")
        equations, guessed_vars = rainbow.rectangular_minrank_attack(
            reduce_dimension=reduce_dimension, verbose=verbose)
        solution = vector([])
        # not implemeted yet
    elif attack_type == 'intersection':
        if verbose:
            print("Mounting the intersection attack...")
        equations, _, _ = rainbow.intersection_attack()
        solution = vector([])
        # not implemeted yet
    if rainbow.q % 2 == 0:
        equations = [delete_powers(eq) for eq in equations]
    return SS, equations, weil_coeff_list, guessed_vars, solution


def get_solution_from_log(log_path, format, N, rainbow=None):
    assert format in ['cms', 'crossbred', 'magma', 'mq', 'wdsat', 'xl']
    with open(log_path, 'r') as f:
        if rainbow != None:
            z = rainbow.F.gens()[0]
            deg = rainbow.ext_deg
            zs = [z ^ i for i in range(deg)]
        variables = []
        found = False
        sol_str = ""
        for line in f.readlines():
            if format == 'cms':
                Nw = N * deg
                if line[0] == 'v':
                    variables += line.strip().split(" ")[1:]
                    if len(variables) < Nw:
                        continue
                    sol = [int((sgn(int(b)) + 1) / 2) for b in variables[:Nw]]
                    parts = [sol[deg * i:deg * i + deg]
                             for i in range(len(sol) / deg)]
                    return vector([linear_combination(bits, zs) for bits in parts])

            if format == 'crossbred':
                if "solution found: " in line:
                    found = True
                    continue
                if found:
                    sol = [ZZ(c) for c in line.strip().strip('][').split(',')]
                    parts = [sol[deg * i:deg * i + deg]
                             for i in range(len(sol) / deg)]
                    return vector([linear_combination(bits, zs) for bits in parts])

            if format == 'magma':
                # find solutions even across multiple lines and take the first one
                if line[:3] == "[ <":
                    found = True
                if found:
                    sol_str += line
                    if ">" in line:
                        end = sol_str.index(">")
                        sol = vector([rainbow.F(x) for x in sol_str[3:end] if x not in [',', ' ', '\n']])
                        parts = [sol[deg * i:deg * i + deg]
                                 for i in range(len(sol) / deg)]
                        return vector([linear_combination(bits, zs) for bits in parts])

            if format == 'mq':
                # possibly includes Weil descent
                if "solution found : " in line:
                    sol = [int(b) for b in line.split(
                        "solution found : ")[1][1:-2].split(", ")]
                    parts = [sol[deg * i:deg * i + deg]
                             for i in range(len(sol) / deg)]
                    return vector([linear_combination(bits, zs) for bits in parts])

            if format == 'wdsat':
                line = line.strip()
                if line == "":
                    continue
                if re.match('^[0-1]*$', line):
                    sol = [ZZ(b) for b in list(line)]
                    parts = [sol[deg * i:deg * i + deg]
                             for i in range(len(sol) / deg)]
                    return vector([linear_combination(bits, zs) for bits in parts])

            if format == 'xl':
                if "  is sol" in line:
                    sol = line.split("  is sol")[0].split(" ")
                    return vector([str_to_elt(rainbow.q, c) for c in sol])
    return None


def check_solution(log_path, solver, solution, N, rainbow, attack_type='differential'):
    if solver in ['cb_orig', 'cb_gpu']:
        log_format = 'crossbred'
    elif solver == 'libfes':
        log_format = 'mq'
    else:
        log_format = solver

    solution_found = get_solution_from_log(
        log_path, format=log_format, N=N, rainbow=rainbow)
    print(f"\n{'First solution found: ' : <25} {solution_found} ")

    if attack_type == 'differential':
        solution_expected = solution[1:-1]
        print(f"{'Expected solution: ' : <25} {solution_expected}\n")
        success = solution_found == solution_expected
        if success:
            print("Attack successful!\n")
        else:
            print("Attack NOT successful. :(\n")
        return success


@ click.command()
@ click.option('--solver', type=click.Choice(['cb_orig', 'cb_gpu', 'cms', 'libfes', 'magma', 'mq', 'wdsat', 'xl'], case_sensitive=False), help='the external solver to be used')
@ click.option('--q', default=16, help='the field order', type=int)
@ click.option('--n', default=0, help='the number of variables', type=int)
@ click.option('--m', default=0, help='the number of equations', type=int)
@ click.option('--o2', default=16, help='the oil subspace dimension', type=int)
@ click.option('--gen_only', default=False, is_flag=True, help='only generate equation systems')
@ click.option('--solve_only', default=False, is_flag=True, help='only solve an existing system and check solutions')
@ click.option('--check_only', default=False, is_flag=True, help='only check solutions')
@ click.option('--inner_hybridation', '-h', default="-1", help='the number of variable that are not guessed in MQ', type=int)
@ click.option('--verbose', '-v', default=False, is_flag=True, help='control the output verbosity')
@ click.option('--reduce_dimension', '-r', default=True, is_flag=True, help='reduce the dimension when possible')
@ click.option('--attack_type', '-t', default='differential', type=click.Choice(['differential', 'minrank', 'intersection'], case_sensitive=False), help='choose attack on Rainbow')
@ click.option('--seed', '-s', default=0, help='the seed for randomness replication', type=int)
@ click.option('--precompiled', default=False, is_flag=True, help='indicates if all relevant solvers are already compiled w.r.t. the parameters')
def main(q, n, m, o2, solver, gen_only, solve_only, check_only, inner_hybridation, verbose, reduce_dimension, attack_type, seed, precompiled):
    if m == 0:
        m = 2*o2
    if n == 0:
        n = 3*o2
    set_random_seed(seed)
    M, N = compute_system_size(q, m, n, o2, attack_type)
    system_folder_path = 'systems'
    Path(system_folder_path).mkdir(parents=True, exist_ok=True)
    log_path = Path("log.txt")
    base_system_name = "rainbow_{}_seed_{}_q_{}_o2_{}_m_{}_n_{}_M_{}_N_{}".format(
        attack_type, seed, q, o2, m, n, M, N)
    cb_gpu_system_path = Path(system_folder_path, base_system_name + '.cb_gpu')
    cb_orig_system_path = Path(system_folder_path, base_system_name + '.cb_orig')
    cnf_system_path = Path(system_folder_path, base_system_name + '.cnf')
    magma_system_path = Path(system_folder_path, base_system_name + '.magma')
    mq_system_path = Path(system_folder_path, base_system_name + '.mq')
    xl_system_path = Path(system_folder_path, base_system_name + '.xl')
    wdsat_system_path = Path(system_folder_path, base_system_name + '.anf')
    setup_path = Path(system_folder_path, base_system_name + '.stp')
    solution_path = Path(system_folder_path, base_system_name + '.sol')
    support = True if attack_type == 'minrank' else False

    if verbose:
        print("Generating Rainbow instance for seed={}, q={}, m={}, n={}, o2={}...".format(
            seed, q, m, n, o2))
    rainbow = Rainbow(q, m, n, o2, support=support, seed=seed)
    if not (solve_only or check_only):
        save_setup(rainbow, setup_path, verbose=verbose)
        SS, equations, weil_coeff_list, guessed_vars, solution = mount_attack(
            rainbow, attack_type, M, N, reduce_dimension=False, verbose=verbose)
        save_solution(solution, solution_path)
        if attack_type == 'differential':
            save_system(file_format='xl', file_path=xl_system_path,
                        rainbow=rainbow, SS=SS, verbose=verbose)
        save_system(file_format='cb_gpu', file_path=cb_gpu_system_path,
                    rainbow=rainbow, SS=SS, equations=equations, weil_coeff_list=weil_coeff_list, verbose=verbose)
        save_system(file_format='cb_orig', file_path=cb_orig_system_path, rainbow=rainbow, equations=equations,
                    guessed_vars=guessed_vars, reduce_dimension=reduce_dimension, verbose=verbose)
        save_system(file_format='cnf', file_path=cnf_system_path, rainbow=rainbow, equations=equations,
                    guessed_vars=guessed_vars, reduce_dimension=reduce_dimension, verbose=verbose)
        save_system(file_format='magma', file_path=magma_system_path, rainbow=rainbow, equations=equations,
                    guessed_vars=guessed_vars, reduce_dimension=reduce_dimension, verbose=verbose)
        save_system(file_format='mq', file_path=mq_system_path, rainbow=rainbow, equations=equations,
                    guessed_vars=guessed_vars, reduce_dimension=reduce_dimension, verbose=verbose)
        save_system(file_format='wdsat', file_path=wdsat_system_path, rainbow=rainbow, equations=equations,
                    guessed_vars=guessed_vars, reduce_dimension=reduce_dimension, verbose=verbose)
    else:
        if verbose:
            print("Skipping the attack equations generation...")
        solution = load_solution(solution_path, q)

    if gen_only or solver == None:
        exit()

    if not check_only:
        equations_path = Path(system_folder_path, base_system_name + f".{get_eq_format(solver)}")
        invoke_solver(solver, equations_path, q, M, N, log_path=log_path,
                      inner_hybridation=inner_hybridation, precompiled=precompiled)

    return check_solution(log_path, solver, solution, N, rainbow)


if __name__ == '__main__':
    main()
