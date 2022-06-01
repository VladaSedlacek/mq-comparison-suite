import click
import pandas as pd
from sage.matrix.symplectic_basis import symplectic_basis_over_field

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


def probability_of_rank(n, m, r, q=2):
    # returns the probability that a random m x n matrix over Fq has rank r
    # computed according to https://www.cs.utoronto.ca/~cvs/coding/random_report.pdf
    rank_r_matrices = gaussian_binomial(n, r)(q=q) * sum([(-1) ^ (r - l) * gaussian_binomial(
        r, l)(q=q) * q ^ (m * l + binomial(r - l, 2)) for l in range(r + 1)])
    all_matrices = q ^ (m * n)
    return rank_r_matrices / all_matrices


def expected_rank(n, m, q=2):
    # returns the expected rank of a random m x n matrix over Fq
    # computed according to https://www.cs.utoronto.ca/~cvs/coding/random_report.pdf
    return sum([r * probability_of_rank(n, m, r, q) for r in range(min(m, n) + 1)])


def expected_optimal_weight(n, m, q=2):
    return sum([expected_rank(i, m, q) for i in range(1, n + 1)]).round()


def max_weight(n, quadratic=True):
    if quadratic:
        return ZZ(n * (n + 1) / 2)
    else:
        return ZZ(n * (n - 1) / 2)


def print_expected_optimal_weights(q=2, n_max=10, m_max=10, n_min=4, m_min=4):
    data = [[f"{expected_optimal_weight(n, m, q)}/{max_weight(n)}"
             for n in range(n_min, n_max + 1)] for m in range(m_min, m_max + 1)]
    df = pd.DataFrame(data)
    df.index = [f"m={m}" for m in range(1, m_max + 1)]
    df.columns = [f"n={n}" for n in range(1, n_max + 1)]
    print(df)


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


def single_symplectic_strategy(MM_bil):
    S_candidates = []
    for M_bil in MM_bil:
        E, S = symplectic_basis_over_field(M_bil)
        assert S * M_bil * S.transpose() == E
        S_candidates.append(S.transpose())
    return S_candidates


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


def get_slices(MM):
    return [Matrix([M[i][i:] for M in MM]) for i in range(MM[0].nrows())]


def guess_optimal_weight(MM):
    return sum([rank(M_slice) for M_slice in get_slices(MM)])


def potential_improvements(MM, i, head=0):
    M_slice = Matrix([M[i][head:] for M in MM])
    zeros = count_zeros_in_vector((Pencil(MM).matrix)[i])
    return corank(M_slice) - zeros


def nonzero_cols(v):
    return [j for j, c in enumerate(v) if c != 0]


def get_vectors_to_try(MM, i, head=0):
    M_slice = Matrix([M[i][head:] for M in MM])
    return [vector([0] * head + list(v)) for v in M_slice.right_kernel().basis()]


def is_in_span(v, vectors_used):
    V = v.parent()
    return v in V.subspace(vectors_used)


def locally_optimal_strategy(MM, extra_tries, quadratic=False, reverse=False, verbose=False):
    if not quadratic:
        reverse = False
    K = MM[0][0, 0].parent()
    n = MM[0].nrows()
    L_total = identity_matrix(K, n)
    for i in range(n):
        head = 0
        if potential_improvements(MM, i) == 0:
            # This might be too restrictive...
            continue
        M_slice = Matrix([M[i][head:] for M in MM])
        if reverse:
            i = n - 1 - i
        if verbose:
            print("\nInvestigating row {}:".format(i))
            print("M_slice:\n{}\nCo-rank: {}".format(
                M_slice, corank(M_slice)))
        vectors_to_try = get_vectors_to_try(MM, i, head)
        vectors_tried = []
        cols_used = []
        extra_counter = 0
        for v in vectors_to_try:
            extra_counter += 1
            if potential_improvements(MM, i) == 0 and extra_counter > extra_tries:
                break
            candidate_cols = nonzero_cols(v)
            if v in vectors_tried or len(candidate_cols) <= 1:
                vectors_tried.append(v)
                continue
            vectors_tried.append(v)
            if verbose:
                print("Trying {} out of possible {}, candidate_cols: {}".format(
                    len(vectors_tried), len(vectors_to_try), candidate_cols))
            for col in candidate_cols:
                if col in cols_used:
                    continue
                L = identity_matrix(K, n)
                L[:, col] = v
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
                    L_total = L_total * L
                    MM = MM_test
                    cols_used.append(col)
                    if verbose:
                        print("Improvement made! New global weight:",
                              global_weight(MM, include_diag=quadratic))
                        if potential_improvements(MM, i) == 0:
                            print("Extra tries actually helped!")
                        print_matrices(MM)
                    break
    print("")
    return L_total


def linear_combination(coefficients, objects):
    assert len(coefficients) == len(objects)
    return sum(c * o for c, o in zip(coefficients, objects))


def thomae_wolf(MM, quadratic=True):
    K = MM[0][0, 0].parent()
    n = MM[0].nrows()
    m = len(MM)
    S = identity_matrix(K, n)
    S[:, 0] = (K ^ n).random_element()
    S[0, 0] = 1

    conditions = []

    for i in range(1, n):
        for M in MM:
            above = list(linear_combination(S.columns()[i - 1], M.rows()))
            below = [linear_combination(S.columns()[i - 1], row)
                     for row in M.rows()]
            conditions.append([a + b for a, b in zip(above, below)])

        C = Matrix(K, conditions)
        if verbose:
            print(f"Solving for column {i}, current rank of C: {C.rank()}")
            print(C.right_kernel(), "\n")
        if C.right_kernel().dimension() == 0:
            break

        for S_col in C.right_kernel():
            assert C * S_col == vector([0] * i * m)

            S_new = copy(S)
            S_new[:, i] = S_col
            if S_new.rank() != n:
                # S[:, i] = C.right_kernel()[-1]
                continue
            else:
                S = S_new

    print(f"\nS (of rank {S.rank()}):")
    print(S, "\n")
    print_matrices(MM, pencilize=True)
    print("")
    print_matrices(bilinear_to_quadratic(
        transform_basis(MM, S)), pencilize=True)

    return S


def combine_matrices(MM):
    K = MM[0][0, 0].parent()
    m = len(MM)
    n = MM[0].nrows()

    MM_new = [linear_combination(v, MM) for v in K ^ m if v != 0]
    MM_new_sorted = sorted(MM_new, key=lambda M: global_weight([M]))
    weights = sorted([global_weight([M]) for M in MM_new])
    print(f"Possible weights: {weights}")
    print(f"Original weights: {[global_weight([M]) for M in MM]}")

    # Chebyshev inequality
    alpha = RR(sqrt(m / 8 * binomial(n, 2))).n(digits=5)
    print(f"Best expected weight: {binomial(n, 2) / 2 - alpha}")

    return MM_new_sorted[:m]


def compare_approaches(MM, quadratic, tries, extra_tries, reverse, verbose=True, width=100):
    K = MM[0][0, 0].parent()
    n = MM[0].nrows()
    print("=" * width + "\n")
    max_weight = ZZ(n * (n - 1) / 2)
    if quadratic:
        max_weight += n
    print("Maximal global weight:", max_weight)
    print("Ranks of slices:", [rank(M_slice) for M_slice in get_slices(MM)])
    print("Estimated optimal weight:", guess_optimal_weight(MM))
    print("With I:")
    if verbose:
        print_matrices(MM)
    print_weight(MM, quadratic=quadratic)

    print("=" * width + "\n")
    print("With locally optimal strategy:")
    L = locally_optimal_strategy(
        MM, quadratic=quadratic, extra_tries=extra_tries, reverse=reverse, verbose=verbose)
    assert L.is_invertible()
    if verbose:
        print_details(MM, L, quadratic=quadratic)
    print_weight(MM, L, quadratic=quadratic)

    print("=" * width + "\n")
    print("With elementary greedy strategy:\n")
    E = elementary_greedy_strategy(MM, tries, quadratic=quadratic)
    if verbose:
        print_details(MM, E, quadratic=quadratic)
    print_weight(MM, E, quadratic=quadratic)

    print("=" * width + "\n")
    print("With Thomae-Wolf strategy:\n")
    S = thomae_wolf(MM, quadratic=quadratic)
    if verbose:
        print_details(MM, S, quadratic=quadratic)
    print_weight(MM, S, quadratic=quadratic)
    print("=" * width)

    return S


@ click.command()
@ click.option('--q', default=2, help='the field order', type=int)
@ click.option('--n', default=10, help='the number of variables', type=int)
@ click.option('--m', default=6, help='the number of equations', type=int)
@ click.option('--o2', default=2, help='the oil subspace dimension', type=int)
@ click.option('-s', '--seed', default=0, help='the seed for randomness replication', type=int)
@ click.option('-t', '--tries', default=100, help='the number of tries for each random strategy', type=int)
@ click.option('-e', '--extra', default=0, help='the maximum number of extra tries for the locally optimal strategy', type=int)
@ click.option('-v', '--verbose', default=False, help='verbosity flag', is_flag=True)
@ click.option('-r', '--reverse', default=False, help='flag for reverse direction in the locally optimal strategy', is_flag=True)
@ click.option('-q', '--quadratic', default=True, help='flag for quadratic strategies instead of bilinear ones', is_flag=True)
@ click.option('-a', '--attack', default=False, help='use the system resulting from the Beullens attack', is_flag=True)
def main(q, n, m, o2, seed, quadratic, tries, extra, reverse, attack, verbose):
    width = 100
    set_random_seed(seed)
    PK, O2, O1, W = Keygen(q, n, m, o2)
    if attack:
        MM = Attack(PK, O2)
    else:
        MM = PK
    MM_bil = [get_polar_form(M) for M in MM]
    print("=" * width)
    print(
        f"The dimensions of the resulting system: n={MM[0].nrows()}, m={len(MM)}")

    # thomae_wolf(MM)
    # combine_matrices(MM)
    # if quadratic:
    #     L = compare_approaches(
    #         MM, quadratic=True, verbose=verbose, tries=tries, extra_tries=extra, reverse=reverse, width=width)
    # else:
    #     L = compare_approaches(MM_bil, quadratic=False,
    #                            verbose=verbose, tries=tries, extra_tries=extra, reverse=reverse, width=width)
    # if verbose:
    #     print("=" * width + "\n")
    #     print("Original quadratic system:")
    #     print_matrices(MM, pencilize=True)
    #     print_weight(MM, quadratic=True)
    #     print("=" * width + "\n")
    #     print("New quadratic system:")
    #     print_details(MM, L, quadratic=True, pencilize=True)
    #     print_weight(MM, L, quadratic=True)


if __name__ == '__main__':
    main()
