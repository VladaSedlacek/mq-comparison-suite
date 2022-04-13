load('rainbow.sage')

# q, n, m, o2 = 2, 12, 8, 4
q, n, m, o2 = 2, 6, 6, 2


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


def find_pair(MM, K):
    found = False
    m = len(MM)
    n = MM[0].ncols()
    Mij = zero_matrix(K, n)
    for i in range(m):
        for j in range(i + 1, m):
            if not found:
                if MM[i].is_invertible() and MM[j].is_invertible():
                    Mij = MM[i] ^ -1 * MM[j]
                    # print(factor(Mij.characteristic_polynomial()))
                    f = poly_sqrt(Mij.characteristic_polynomial())
                    if f.is_irreducible():
                        found = True
                        return Mij, i, j


set_random_seed(2)
K = GF(q)
I = identity_matrix(K, n)
PK, O2, O1, W = Keygen(q, n, m, o2)
MM = [get_polar_form(M) for M in PK]

Mij, i, j = find_pair(MM, K)

fact = factor(Mij.characteristic_polynomial())
f = poly_sqrt(Mij.characteristic_polynomial())
print("char poly:", fact, "\nsqrt:     ", factor(f))
Cf = companion_matrix(f, format='bottom')
# C = block_diagonal_matrix(Cf, Cf.transpose())
C = block_diagonal_matrix(Cf, Cf)
print(C, "\n")
assert C.is_similar(block_diagonal_matrix(Cf, Cf))
b, P = C.is_similar(Mij, transformation=True)
assert b
assert P.is_invertible()

Mi = (P.transpose() * MM[i] * P)[n / 2:, :n / 2]
Mj = (P.transpose() * MM[j] * P)[n / 2:, :n / 2]
MiD = block_matrix([[0, Mi], [Mi, 0]])
MjD = block_matrix([[0, Mj], [Mj, 0]])
assert P.transpose() * MM[i] * P == MiD
assert P.transpose() * MM[j] * P == MjD
assert C.is_similar(block_diagonal_matrix(Mi ^ -1 * Mj, Mj * Mi ^ -1))
print(MiD, "\n")
print(MjD, "\n")

print("Before:")
print_matrices(MM)
print("Total weight:", total_weight(MM, I))
print("Global weight:", global_weight(MM, I))

print("\nWith P:")
print_matrices(transform_basis(MM, P))
print("Total weight:", total_weight(MM, P))
print("Global weight:", global_weight(MM, P))

print("\nAt random:")
record = total_weight(MM, P)
for i in range(100):
    U = random_matrix(K, n)
    if not U.is_invertible():
        continue
    if total_weight(MM, U) < record:
        record = total_weight(MM, U)
        print("Transformation matrix:\n", U, "\n")
        print_matrices(transform_basis(MM, U))
        print("Total weight:", total_weight(MM, U))
        print("Global weight:", global_weight(MM, U))
