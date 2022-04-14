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


def transform_basis(PK, U):
    return [U.transpose() * M * U for M in PK]


def total_weight(PK, U):
    return sum([hamming_weight(M) for M in transform_basis(PK, U)])


def global_weight(PK, U):
    weight = 0
    for i in range(PK[0].nrows()):
        for j in range(i + 1, PK[0].ncols()):
            all_zeros = True
            for M in transform_basis(PK, U):
                if M[i, j] != 0:
                    all_zeros = False
                    break
            if not all_zeros:
                weight += 1
    return weight


def print_matrices(PK):
    for i in range(PK[0].nrows()):
        for M in PK:
            print(list(M[i]), end="\t")
        print("")


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


def print_weights(MM, U, show_total_weight=False):
    if not show_total_weight:
        print("Global weight: {}".format(global_weight(MM, U)))
    else:
        print("Global weight: {}, total weight: {}".format(
            global_weight(MM, U), total_weight(MM, U)))


def compare_approaches(MM, tries=100, show_matrices=True, show_total_weight=False):
    K = MM[0][0, 0].parent()
    I = identity_matrix(K, n)

    print("With I:")
    if show_matrices:
        print_matrices(MM)
    print_weights(MM, I, show_total_weight)

    print("\nWith symplectic basis for two matrices:")
    S = find_symplectic_for_two(MM)
    print("Transformation matrix:\n{}\n".format(S))
    if show_matrices:
        print_matrices(transform_basis(MM, S))
    print_weights(MM, S, show_total_weight)

    print("\nAt random:")
    R = find_best_random(MM, tries)
    print("Transformation matrix:\n{}\n".format(R))
    if show_matrices:
        print_matrices(transform_basis(MM, R))
    print_weights(MM, R, show_total_weight)


# q, n, m, o2 = 2, 12, 8, 4
q, n, m, o2 = 2, 6, 6, 2
set_random_seed(2)
PK, O2, O1, W = Keygen(q, n, m, o2)
MM = [get_polar_form(M) for M in PK]
compare_approaches(MM, tries=1000)
