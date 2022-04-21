import click
load('rainbow.sage')


def Attack(PK, O2=None):
    m = len(PK)
    n = PK[0].ncols()
    K = PK[0][0, 0].parent()
    q = K.order()
    basis_Fn = (K**n).basis()

    global attempts
    attempts = 0

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


def get_polar_form(Q):
    return Q + Q.transpose()


def bilinear_to_quadratic(MM):
    for M in MM:
        Make_UD(M)
    return MM


def transform_basis(MM, U):
    return [U.transpose() * M * U for M in MM]


def global_weight(MM, U=None, max_row=None, max_col=None, include_diag=False):
    if U is not None:
        MM = transform_basis(MM, U)
    if include_diag:
        MM = bilinear_to_quadratic(MM)
    if max_row is None:
        max_row = MM[0].nrows()
    if max_col is None:
        max_col = MM[0].ncols()
    weight = 0
    for i in range(max_row):
        min_col = i + 1 - int(include_diag)
        for j in range(min_col, max_col):
            if Pencil(MM).matrix[i, j] != 0:
                weight += 1
    return weight


def print_weight(MM, U=None, quadratic=False):
    print("Global weight: {}".format(
        global_weight(MM, U, include_diag=quadratic)))


def print_matrices(MM, pencilize=True):
    if pencilize:
        print(Pencil(MM).matrix)
    else:
        for i in range(MM[0].nrows()):
            for M in MM:
                print(list(M[i]), end="\t")
            print("")


def print_details(MM, U, quadratic=False, pencilize=True):
    print("Transformation matrix:\n{}\n".format(U))
    transformed = transform_basis(MM, U)
    if quadratic:
        transformed = bilinear_to_quadratic(transformed)
    print_matrices(transformed, pencilize)


def poly_sqrt(poly):
    lst = list(factor(poly))
    result = 1
    for f, e in lst:
        assert e % 2 == 0
        result *= f ^ (ZZ(e / 2))
    return result


def find_pair(MM):
    m = len(MM)
    n = MM[0].ncols()
    K = MM[0][0, 0].parent()
    Mij = zero_matrix(K, n)
    for i in range(m):
        for j in range(i + 1, m):
            if MM[i].is_invertible() and MM[j].is_invertible():
                Mij = MM[i].inverse() * MM[j]
                f = poly_sqrt(Mij.characteristic_polynomial())
                if f.is_irreducible():
                    return Mij, i, j
    print("No good ratio of matrices found!")
    return None, None, None


def find_best_random(MM, tries=100):
    n = MM[0].ncols()
    K = MM[0][0, 0].parent()
    record_U = identity_matrix(K, n)
    record_weight = global_weight(MM)
    for _ in range(tries):
        U = random_matrix(K, n)
        if not U.is_invertible():
            continue
        if global_weight(MM, U) < record_weight:
            record_weight = global_weight(MM, U)
            record_U = U
    return record_U


def find_symplectic_for_two(MM, verbose=False, checks=False):
    # TODO: generalize this to linear combinations
    Mij, i, j = find_pair(MM)
    if Mij == None:
        K = MM[0][0, 0].parent()
        n = MM[0].nrows()
        return identity_matrix(K, n)
    fact = factor(Mij.characteristic_polynomial())
    f = poly_sqrt(Mij.characteristic_polynomial())
    if verbose:
        print("char poly:", fact, "\nsqrt:     ", factor(f))
    Cf = companion_matrix(f, format='bottom')
    # C = block_diagonal_matrix(Cf, Cf.transpose())
    C = block_diagonal_matrix(Cf, Cf)
    if verbose:
        print(C, "\n")
    b, U = C.is_similar(Mij, transformation=True)
    if checks:
        assert C.is_similar(block_diagonal_matrix(Cf, Cf))
        assert b
        assert U.is_invertible()
        Mi = (U.transpose() * MM[i] * U)[n / 2:, :n / 2]
        Mj = (U.transpose() * MM[j] * U)[n / 2:, :n / 2]
        MiD = block_matrix([[0, Mi], [Mi, 0]])
        MjD = block_matrix([[0, Mj], [Mj, 0]])
        assert U.transpose() * MM[i] * U == MiD
        assert U.transpose() * MM[j] * U == MjD
        assert C.is_similar(block_diagonal_matrix(Mi ^ -1 * Mj, Mj * Mi ^ -1))
    if verbose:
        print(MiD, "\n")
        print(MjD, "\n")
    return U


def elementary_improvement(MM, i, j, s=1, quadratic=False, verbose=False):
    # try adding the j-th row scaled by s to the i-th row
    K = MM[0][0, 0].parent()
    n = MM[0].nrows()
    E = elementary_matrix(K, n, row1=i, row2=j, scale=s)
    old_weight = global_weight(MM, include_diag=quadratic)
    new_weight = global_weight(MM, E, include_diag=quadratic)
    improved = old_weight > new_weight
    if verbose and improved:
        print("weights:", old_weight, new_weight)
    return improved, E


def elementary_greedy_strategy(MM, tries=100, quadratic=False, I_start=None):
    K = MM[0][0, 0].parent()
    n = MM[0].nrows()
    E_total = identity_matrix(K, n)
    if I_start is not None:
        MM = [I_start.transpose() * M * I_start for M in MM]
        E_total = I_start
    for _ in range(tries):
        i = Integers(n).random_element()
        j = Integers(n).random_element()
        if i == j:
            continue
        improved, E = elementary_improvement(MM, i, j)
        if improved:
            E_total = E_total * E
            MM = [E.transpose() * M * E for M in MM]
            MM = bilinear_to_quadratic(MM)
    return E_total


def elementary_greedy_strategy_iterated(MM, tries=100, starts=5):
    n = MM[0].ncols()
    K = MM[0][0, 0].parent()
    record_U = identity_matrix(K, n)
    record_weight = global_weight(MM)
    for _ in range(starts):
        from sage.matrix.constructor import random_unimodular_matrix
        I_start = random_unimodular_matrix(
            sage.matrix.matrix_space.MatrixSpace(K, n))
        if not I_start.is_invertible():
            print("nope")
            continue
        U = elementary_greedy_strategy(MM, tries=tries, I_start=I_start)
        if global_weight(MM, U) < record_weight:
            record_weight = global_weight(MM, U)
            record_U = U
    return record_U


class Pencil(object):
    """a class representing a pencil of bilinear forms"""

    def __init__(self, MM, tries=100):
        self.MM = MM
        self.m = len(MM)
        self.n = MM[0].ncols()
        self.K = MM[0][0, 0].parent()
        self.I = identity_matrix(self.K, self.n)
        self.PR = PolynomialRing(self.K, 'z')
        self.z = self.PR.gen()
        self.module = FreeModule(self.K, len(MM))
        self.matrix = sum([M * self.z ^ i for i, M in enumerate(MM)])
        self.tries = tries

    def custom_reduce(self, restart=True):
        MMP = self.matrix
        E_total = self.I
        E_best = E_total
        tries = self.tries
        while tries > 0:
            if restart:
                E_total = self.I
            perm = Permutations(self.n).random_element()
            for col in perm.action(MMP.columns()):
                for deg in range(self.m):
                    perm = Permutations(self.n).random_element()
                    pivot_poly = 0
                    pivot_index = -1
                    for i, poly in perm.action(list(enumerate(col))):
                        if get_coeff(poly, deg) != 0:
                            tries = tries - 1
                            if pivot_poly == 0:
                                pivot_poly = poly
                                pivot_index = i
                            else:
                                poly += 1 * pivot_poly
                                # add the j-th row scaled by 1 to the i-th row
                                E = elementary_matrix(
                                    self.K, self.n, row1=i, row2=pivot_index, scale=1)
                                E_total = E_total * E
            if global_weight(self.MM, E_best) > global_weight(self.MM, E_total):
                E_best = E_total
        return E_best

    def reduce(self, deg):
        _, U = self.matrix.reduced_form(transformation=True)
        Uc = coefficient_matrix(U, deg)
        return self.invertible_or_identity(Uc)

    def hermite(self, deg):
        _, U = self.matrix.hermite_form(transformation=True)
        Uc = coefficient_matrix(U, deg)
        return self.invertible_or_identity(Uc)

    def smith_right(self, deg):
        D, U, V = self.matrix.smith_form()
        Vc = coefficient_matrix(U, deg)
        return self.invertible_or_identity(Vc)

    def smith_left(self, deg):
        D, U, V = self.matrix.smith_form()
        Uc = coefficient_matrix(U, deg)
        return self.invertible_or_identity(Uc)

    def invertible_or_identity(self, U):
        if U.is_invertible():
            return U
        return self.I


def all_coeffs(poly, max_deg):
    return list(poly) + (max_deg + 1 - len(list(poly))) * [0]


def get_coeff(poly, deg):
    return all_coeffs(poly, deg)[deg]


def coefficient_matrix(poly_matrix, deg):
    rows = [vector(get_coeff(poly, deg) for poly in row)
            for row in poly_matrix]
    return Matrix(rows)


def corank(M):
    return M.ncols() - M.rank()


def count_zeros_in_vector(v):
    return sum([int(vi == 0) for vi in v])


def locally_optimal_strategy(MM, quadratic=False, reverse=False, try_all=True, verbose=False):
    if not quadratic:
        reverse = False
    K = MM[0][0, 0].parent()
    n = MM[0].nrows()
    L_total = identity_matrix(K, n)
    for i in range(n):
        if reverse:
            i = n - 1 - i
        M_across = Matrix([M[i] for M in MM])
        if verbose:
            print("\nInvestigating row {}:".format(i))
            print("M_across:\n{}\nCo-rank: {}".format(
                M_across, corank(M_across)))
        ker = M_across.right_kernel()
        if try_all:
            vectors_to_try = ker
        else:
            vectors_to_try = ker.basis_matrix()
        zeros_in_row = count_zeros_in_vector((Pencil(MM).matrix)[i])
        potential_improvements = ker.dimension() - zeros_in_row
        while potential_improvements > 0:
            if verbose:
                print("Potential improvements:", potential_improvements)
            for ker_vec in vectors_to_try:
                candidate_cols = [j for j, c in enumerate(ker_vec) if c != 0]
                if len(candidate_cols) == 1:
                    potential_improvements -= 1
                    continue
                if verbose:
                    print("Candidate cols:", candidate_cols)
                for col in candidate_cols:
                    L = identity_matrix(K, n)
                    L[:, col] = ker_vec
                    if quadratic:
                        MM_test = bilinear_to_quadratic(
                            transform_basis(MM, L))
                        condition = global_weight(MM_test, include_diag=True) < global_weight(
                            MM, include_diag=True)
                    else:
                        MM_test = transform_basis(MM, L)
                        condition = global_weight(
                            MM, L, max_row=i) < global_weight(MM, max_row=i)
                    if condition:
                        assert L.is_invertible()
                        L_total = L_total * L
                        MM = MM_test
                        if verbose:
                            print("Improvement made! New global weight:",
                                  global_weight(MM, include_diag=quadratic))
                            print_matrices(MM)
                        break
                potential_improvements -= 1
        if verbose:
            print("No potential improvements, moving on...")
    print("")
    return L_total


def compare_approaches(MM, tries=100, quadratic=False, verbose=True, width=100):
    K = MM[0][0, 0].parent()
    n = MM[0].nrows()
    print("=" * width + "\n")
    max_weight = ZZ(n * (n - 1) / 2)
    if quadratic:
        max_weight += n
    print("Maximal global weight:", max_weight)
    print("With I:")
    if verbose:
        print_matrices(MM)
    print_weight(MM, quadratic=quadratic)

    print("=" * width + "\n")
    print("With locally optimal strategy:")
    L = locally_optimal_strategy(
        MM, quadratic=quadratic, reverse=False, verbose=False)
    if verbose:
        print_details(MM, L, quadratic=quadratic)
    print_weight(MM, L, quadratic=quadratic)

    print("=" * width + "\n")
    print("With elementary greedy strategy:\n")
    E = elementary_greedy_strategy(MM, tries, quadratic=quadratic)
    if verbose:
        print_details(MM, E, quadratic=quadratic)
    print_weight(MM, E, quadratic=quadratic)
    print("=" * width)

    return L


@ click.command()
@ click.option('--q', default=2, help='the field order', type=int)
@ click.option('--n', default=10, help='the number of variables', type=int)
@ click.option('--m', default=6, help='the number of equations', type=int)
@ click.option('--o2', default=2, help='the oil subspace dimension', type=int)
@ click.option('-s', '--seed', default=0, help='the seed for randomness replication', type=int)
@ click.option('-t', '--tries', default=100, help='the number of tries for each random strategy', type=int)
@ click.option('-v', '--verbose', default=False, help='verbosity flag', is_flag=True)
@ click.option('-q', '--quadratic', default=False, help='flag for quadratic strategies instead of bilinear ones', is_flag=True)
def main(q, n, m, o2, seed, tries, verbose, quadratic):
    width = 100
    set_random_seed(seed)
    PK, O2, O1, W = Keygen(q, n, m, o2)
    MM = [get_polar_form(M) for M in PK]
    SS = Attack(PK, O2)
    SS_bil = [get_polar_form(S) for S in SS]
    print("=" * width)

    if quadratic:
        L = compare_approaches(
            SS, quadratic=True, verbose=verbose, tries=tries, width=width)
    else:
        L = compare_approaches(SS_bil, quadratic=False,
                               verbose=verbose, tries=tries, width=width)
    if verbose:
        print("=" * width + "\n")
        print("Original quadratic system:")
        print_matrices(SS, pencilize=True)
        print_weight(SS, quadratic=True)
        print("=" * width + "\n")
        print("New quadratic system:")
        print_details(SS, L, quadratic=True, pencilize=True)
        print_weight(SS, L, quadratic=True)


if __name__ == '__main__':
    main()
