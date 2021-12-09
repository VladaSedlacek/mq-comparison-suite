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
        self.k = find_max_k(self.m, self.o2, verbose=False)
        self.reduced = self.q % 2 == 0 and self.n % 2 == 1
        self.V = VectorSpace(F, n)
        self.R = PolynomialRing(F, ['x%s' % p for p in range(
            1, n + 1)], order="neglex")
        self.R.inject_variables()
        self.xx = vector(self.R.gens()[:n])
        self.FF = self.construct_central_map()
        self.T, self.S, self.PP, self.MM = self.hide_central_map()

    def construct_central_map(self):
        m, n = self.m, self.n
        FF = []
        for _ in range(o2):
            Q = Matrix(self.F, n)
            for i in range(n):
                for j in range(n):
                    if i >= n - o2 and j >= n - m:
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
            if T.rank() == m:
                break
        PP = [S * T.transpose() * Q * T for Q in self.FF]
        MM = [get_polar_form(P) for P in PP]
        assert MM == [get_polar_form(P) for P in PP]
        return T, S, PP, MM

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
