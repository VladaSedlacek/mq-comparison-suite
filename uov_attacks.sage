from itertools import product

def construct_central_map(F, m, n):
    FF = []
    for _ in range(m):
        Q = Matrix(F, n)
        for i in range(n):
            for j in range(n):
                if i >= n-m and j >= n-m:
                    continue
                Q[i,j] = F.random_element()
        FF.append(Q)
    return FF

def get_polar_form(Q):
    return Q+Q.transpose()

def setup(q, m, n, variables=True, verbose=False):
    assert m < n
    F = GF(q)
    V = VectorSpace(F, n)

    if variables:
        R = PolynomialRing(F, ['x%s'%p for p in range(1,n+1)] + ['y%s'%p for p in range(1,n+1)])
        R.inject_variables()
        xx = vector(R.gens()[:n])
        yy = vector(R.gens()[n:])
        W = VectorSpace(F, m)

    FF = construct_central_map(F, m, n)
    T = random_matrix(F, n)
    assert T.rank() == n
    PP = [T.transpose()*Q*T for Q in FF]
    MM = [T.transpose()*get_polar_form(Q)*T for Q in FF]
    assert MM == [get_polar_form(P) for P in PP]

    O2_basis = [V([0]*n) for _ in range(m)]
    for i in range(m):
        O2_basis[i][i+n-m] = 1
    for o in O2_basis:
        for Q in FF:
            assert o*Q*o ==0
    O_basis = [T.inverse()*o for o in O2_basis]
    for o in O_basis:
        for P in PP:
            assert o*P*o==0
    for i in range(m):
        for j in range(i, m):
            for M in MM:
                assert O_basis[i]*M*O_basis[j]==0
    O = V.subspace(O_basis)
    return xx, FF, PP, MM, T, O

def kipnis_shamir(m, n, MM, verbose=False):
    assert n == 2 * m
    inv_subspaces = []
    for i in range(m):
        for j in range(i+1, m):
            try:
                Mij = MM[i] * MM[j].inverse()
            except ZeroDivisionError:
                continue
            if verbose:
                print("i,j:", i,j)
                print("Mi,Mj:", MM[i], MM[j])
                print("Mij:", Mij)
            poly = Mij.characteristic_polynomial().radical()
            if poly(Mij) == 0:
                continue
            inv_subspaces.append(poly(Mij).kernel())
    return inv_subspaces

def intersection_attack(m, n, MM, verbose=False):
    assert n < 3 * m
    found = 0
    equations = []
    for i in range(m):
        for j in range(i+1, m):
            try:
                Mi_inv = MM[i].inverse()
                Mj_inv = MM[j].inverse()
                found = 1
            except ZeroDivisionError:
                continue
            if verbose:
                print("i,j:", i,j)
                print("Mi_inv, Mj_inv:\n", Mi_inv, "\n\n", Mj_inv, Mj_inv*xx)
            if found == 1:
               u = Mi_inv*xx
               v = Mj_inv*xx
               for P in PP:
                   equations.append(u * P * u)
                   equations.append(v * P * v)
               for M in MM:
                   equations.append(u * M * v)
               return (equations, i, j)
            return ([], i, j)


def intersection_attack_advanced(m, n, MM, verbose=False):
    assert n > 2*m and n < 3*m
    k = 2
    while True:
        if verbose:
            print(k, RR(n/m), RR((2*k-1)/(k-1)))
        if n >= (2*k-1)/(k-1) * m:
            k -= 1
            break
        k += 1
    if verbose:
        print("k:", k)
    LL = []
    for i in range(k):
        L = Matrix(GF(q),n)
        while not L.is_invertible():
            for j in range(len(MM)):
                L += GF(q).random_element() * MM[j]
        LL.append(L)
    if verbose:
        print(LL)
    equations = []
    for i in range(k):
        u = LL[i].inverse() * xx
        for P in PP:
            equations.append(u * P * u)
        for j in range(i+1, k): 
            v = LL[j].inverse() * xx
            for M in MM:
                equations.append(u * M * v)
    return (equations, None, None)

q = 11
m = 3
n = 7
xx, FF, PP, MM, T, O = setup(q, m, n)

# print(kipnis_shamir(m,n,MM))
print(intersection_attack(m,n,MM))
print(intersection_attack_advanced(m,n,MM))