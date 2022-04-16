import click
load('rainbow.sage')


def get_polar_form(Q):
    return Q + Q.transpose()


def hamming_weight(M):
    weight = 0
    for i in range(M.nrows()):
        for j in range(i + 1, M.ncols()):
            if M[i, j] != 0:
                weight += 1
    return weight


def transform_basis(MM, U):
    return [U.transpose() * M * U for M in MM]


def total_weight(MM, U):
    return sum([hamming_weight(M) for M in transform_basis(MM, U)])


def global_weight(MM, U):
    weight = 0
    for i in range(MM[0].nrows()):
        for j in range(i + 1, MM[0].ncols()):
            all_zeros = True
            for M in transform_basis(MM, U):
                if M[i, j] != 0:
                    all_zeros = False
                    break
            if not all_zeros:
                weight += 1
    return weight


def print_matrices(MM, pencilize=True):
    if pencilize:
        print(Pencil(MM).matrix)
    else:
        for i in range(MM[0].nrows()):
            for M in MM:
                print(list(M[i]), end="\t")
            print("")


def print_details(MM, U, pencilize=True):
    print("Transformation matrix:\n{}\n".format(U))
    print_matrices(transform_basis(MM, U), pencilize)


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
    I = identity_matrix(K, n)
    record_U = I
    record_weight = global_weight(MM, I)
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


def elementary_improvement(MM, i, j, s=1, verbose=False):
    # try adding the j-th row scaled by s to the i-th row
    K = MM[0][0, 0].parent()
    n = MM[0].nrows()
    I = identity_matrix(K, n)
    E = elementary_matrix(K, n, row1=i, row2=j, scale=s)
    old_weight = global_weight(MM, I)
    new_weight = global_weight(MM, E)
    improved = old_weight > new_weight
    if verbose and improved:
        print("weights:", old_weight, new_weight)
    return improved, E


def elementary_greedy_strategy(MM, tries=100, I_start=None):
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
    return E_total


def elementary_greedy_strategy_iterated(MM, tries=100, starts=5):
    n = MM[0].ncols()
    K = MM[0][0, 0].parent()
    I = identity_matrix(K, n)
    record_U = I
    record_weight = global_weight(MM, I)
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


def print_weights(MM, U, show_total_weight=False):
    if not show_total_weight:
        print("Global weight: {}".format(global_weight(MM, U)))
    else:
        print("Global weight: {}, total weight: {}".format(
            global_weight(MM, U), total_weight(MM, U)))


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


def compare_approaches(MM, tries=100, verbose=True, show_total_weight=False):
    K = MM[0][0, 0].parent()
    n = MM[0].nrows()
    I = identity_matrix(K, n)
    print("Maximal global weight:", ZZ(n * (n - 1) / 2))
    print("With I:")
    if verbose:
        print_matrices(MM)
    print_weights(MM, I, show_total_weight)

    print("\nWith symplectic basis for two matrices:")
    S = find_symplectic_for_two(MM)
    if verbose:
        print_details(MM, S)
    print_weights(MM, S, show_total_weight)

    print("\nAt random:")
    R = find_best_random(MM, tries)
    if verbose:
        print_details(MM, R)
    print_weights(MM, R, show_total_weight)

    print("\nWith elementary greedy strategy:")
    E = elementary_greedy_strategy(MM, tries)
    if verbose:
        print_details(MM, E)
    print_weights(MM, E, show_total_weight)

    print("\nWith custom pencil strategy:")
    pencil = Pencil(MM, tries=tries)
    E = pencil.custom_reduce(restart=True)
    if verbose:
        print_details(MM, E)
    print_weights(MM, E, show_total_weight)


@ click.command()
@ click.option('--q', default=2, help='the field order', type=int)
@ click.option('--n', default=10, help='the number of variables', type=int)
@ click.option('--m', default=6, help='the number of equations', type=int)
@ click.option('--o2', default=2, help='the oil subspace dimension', type=int)
@ click.option('-s', '--seed', default=0, help='the seed for randomness replication', type=int)
@ click.option('-t', '--tries', default=100, help='the number of tries for each random strategy', type=int)
@ click.option('-v', '--verbose', default=False, help='the number of tries for each random strategy', is_flag=True)
def main(q, n, m, o2, seed, tries, verbose):
    set_random_seed(seed)
    PK, O2, O1, W = Keygen(q, n, m, o2)
    MM = [get_polar_form(M) for M in PK]
    compare_approaches(MM, verbose=verbose, tries=tries)


if __name__ == '__main__':
    main()
