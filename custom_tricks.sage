import itertools
load('rainbow.sage')


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


def print_matrices(MM, pencilize=True):
    if pencilize:
        print(Pencil(MM).matrix)
    else:
        for i in range(MM[0].nrows()):
            for M in MM:
                print(list(M[i]), end="\t")
            print("")


def monomials(vars, degree):
    '''Compute all monomials of a certain degree'''
    if degree < 0:
        return
    for comb in itertools.combinations_with_replacement(vars, degree):
        u = 1
        for var in comb:
            u *= var
        yield u
    return


def delete_powers(eq):
    return sum([radical(mon) for mon in eq.monomials()])


def macaulay(eqs, vars, D=2):
    mons = monomials(vars, D - 2)
    for f in eqs:
        for u in mons:
            uf = u * f
            if uf != 0 and uf not in eqs:
                eqs.append(uf)
    s = Sequence([delete_powers(eq) for eq in eqs])
    Mac, cols = s.coefficient_matrix()
    return Mac, cols


def function_of_mon_sum(n, m, l, r=2):
    assert l <= m
    print(f"Mac_2 has {l} rows and {RR(binomial(n,2) - r).n(digits=5)} cols")
    print(
        f"Mac_3 has {n*l} rows and {RR(binomial(n,2) - r + binomial(n,3) - r).n(digits=5)} cols")
    return RR(binomial(binomial(n, 2), r) * binomial(m, l) * 2 ^ (-l * (r - 1))).n(digits=5)


def is_sublist(ls1, ls2):
    def get_all_in(one, another):
        for element in one:
            if element in another:
                yield element

    for x1, x2 in zip(get_all_in(ls1, ls2), get_all_in(ls2, ls1)):
        if x1 != x2:
            return False

    return True


def find_coherent_mac_rows(Mac, l, r=2):
    assert r == 2
    viable_rows_list = []
    # for col1, col2 in itertools.combinations(Mac.columns(), r):
    for i in range(Mac.ncols()):
        for j in range(i + 1, Mac.ncols()):
            col1 = Mac.columns()[i]
            col2 = Mac.columns()[j]
            diff = [c1 == c2 for c1, c2 in zip(col1, col2)]
            score = sum(diff)
            if score >= l:
                # print(i, j)
                viable_rows = ([i for i, b in enumerate(diff) if b == True])
                viable_rows_list.append(viable_rows)
    viable_rows_list = sorted(viable_rows_list, key=lambda x: -len(x))
    intersect = viable_rows_list[0]
    cut = 1
    while len(intersect) > l and cut < len(viable_rows_list):
        intersect = sorted(list(set.intersection(
            *map(set, viable_rows_list[:cut]))))
        cut += 1
    # print(viable_rows_list)
    print(f"\n{intersect}\n")
    MacR = Mac[intersect]
    # print(MacR)

    pairs = []
    for i in range(MacR.ncols()):
        for j in range(i + 1, MacR.ncols()):
            col1 = MacR.columns()[i]
            col2 = MacR.columns()[j]
            diff = [c1 == c2 for c1, c2 in zip(col1, col2)]
            score = sum(diff)
            if score >= l:
                pairs.append((i, j))
                # viable_rows = ([i for i, b in enumerate(diff) if b == True])
                # if is_sublist(viable_rows, intersect):
                #     print(i, j)
    return intersect, MacR, pairs


# for n in range(30, 40):
#     print(n, function_of_mon_sum(n, n, n - 10))
n = 8
m = 5
l = 8
r = 2
PR = PolynomialRing(GF(2), n, 'x')
PR.inject_variables(verbose=False)
xx = vector(PR, [PR.gens()])
set_random_seed(0)

MM = [random_matrix(GF(2), n) for _ in range(m)]
for M in MM:
    Make_UD(M)
eqs = [xx * M * xx for M in MM]

# print_matrices(MM, pencilize=False)
# for eq in eqs:
#     print(eq)


# Mac, cols = macaulay(eqs, xx, 2)

# print("Macaulay matrix:\n", cols.transpose(), "\n")
# print(Mac)
# print(Mac.dimensions())

# intersect, MacR, pairs = find_coherent_mac_rows(Mac, l)

# print(MacR, "\n")
# new_cols = []
# to_delete = [j for i, j in pairs]
# MacR = MacR.delete_columns(to_delete)

# # print(cols.transpose())
# print(MacR)
# print(MacR.dimensions())

def count_frequency(eq, n):
    frequency = {"0": 0, "1": 0}
    for v in GF(2) ^ n:
        frequency[str(eq(*v))] += 1
    return frequency


def split_eq(eq):
    eq_quad = sum(
        [mon for coeff, mon in delete_powers(eq) if mon.degree() == 2])
    eq_lin = delete_powers(eq) - eq_quad
    return eq_quad, eq_lin


def count_partial_frequencies(eq, n):
    eq_quad, eq_lin = split_eq(eq)
    print(eq_quad, "\t", eq_lin)
    print(count_frequency(eq, n))
    print(count_frequency(eq_quad, n))
    print(count_frequency(eq_lin, n))


def count_nondeg_frequencies(n):
    if n % 2 == 1:
        N = (n - 1) / 2
        return sum([binomial(N, 2 * r) * 3 ^ (N - 2 * r) for r in [0..floor(N / 2)]])
    else:
        return count_nondeg_frequencies(n - 1)


def possible_frequencies(n):
    if n % 2 == 1:
        possible = [1 / 2]
        for l in range(0, (n - 1) / 2 + 1):
            freq = count_nondeg_frequencies(2 * l + 1) / 2 ^ (2 * l)
            print(l, freq)
            possible.append(freq)
        return possible


# eqs = [xx * M * xx for M in MM]
# for eq in eqs:
    # count_partial_frequencies(eq, n)

while True:
    S = random_matrix(GF(2), n)
    if S.is_invertible():
        break

print("=" * 78)

# eqs = [xx * S.transpose() * M * S * xx for M in MM]
# for eq in eqs:
# count_partial_frequencies(eq, n)

# print(possible_frequencies(11))


MM = [random_matrix(GF(2), n) for _ in range(m)]
for M in MM:
    Make_UD(M)
    for i in range(m):
        for j in range(i, m):
            M[i, j] = 0

for M in MM:
    print(M * S, "\n")
