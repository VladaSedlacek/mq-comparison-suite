#!/usr/bin/env sage

from itertools import product
from pathlib import Path
import click
import subprocess as sp
from invoke_solver import invoke_solver
from config_utils import declare_paths, get_eq_path, get_sol_path, use_weil
load("equation_utils.sage")


# multivariate util functions

def get_polar_form(Q):
    return Q + Q.transpose()


def Eval(F, x):
    # evaluate Multivariate map and Differential
    return vector([x * M * x for M in F])


def Differential(F, x, y):
    return vector([(x * M * y) + (y * M * x) for M in F])


def Make_UD(M):
    # make a matrix upper diagonal
    n = M.ncols()
    K = M[0, 0].parent()
    for i in range(n):
        for j in range(i + 1, n):
            M[i, j] += M[j, i]
            M[j, i] = K(0)


class Rainbow():
    """A class for the Rainbow scheme."""

    def __init__(self, q, m, n, o2, seed=0, debug=True):
        assert 2 < o2 and o2 < m and m < n
        self.seed = seed
        self.debug = debug
        self.q = q
        F = GF(q)
        self.ext_deg = F.modulus().degree()
        self.F = F
        self.m = m
        self.n = n
        self.o2 = o2
        self.V = VectorSpace(F, n)
        self.V2 = VectorSpace(F, m)
        self.order = "degrevlex"
        self.R = PolynomialRing(F, ['x%s' % p for p in range(1, n + 1)], order=self.order)
        self.xx = vector(self.R.gens())
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

            if alpha == 0:
                NewSol = xt
            else:
                NewSol = xt + 1 / alpha.sqrt() * YSol
            if NewSol[-1] != 0:
                NewSol = NewSol / NewSol[-1]
            # assert NewSol == Sol
            # assert Eval(SS_orig, NewSol) == vector(m * [0])

            # handle indexing from 0/1
            affine_vars = vector(list(self.xx[: n-m-2]) + [1])
            equations = [affine_vars * s * affine_vars for s in SS]
        # TODO odd char
        return equations, Sol


def save_setup(rainbow, setup_path, verbose=False):
    if setup_path.is_file() and verbose:
        print("The file {} already exists!".format(str(setup_path)))
        return
    with open(setup_path, 'w') as f:
        f.write("# seed:" + str(rainbow.seed) + "\n")
        f.write("# O1:" + str(rainbow.O1) + "\n\n")
        f.write("# O2:" + str(rainbow.O2) + "\n\n")
        f.write("# W:" + str(rainbow.W) + "\n\n")


def compute_system_size(q, m, n, attack_type='differential'):
    '''Return the number of equations and the number of variables'''
    if attack_type == 'differential':
        if q % 2 == 0:
            return m - 1, n - m - 2
        return m, n - m - 1


@ click.command()
@ click.option('--q', default=16, help='the field order', type=int)
@ click.option('--n', default=0, help='the number of variables', type=int)
@ click.option('--m', default=0, help='the number of equations', type=int)
@ click.option('--o2', default=16, help='the oil subspace dimension', type=int)
@ click.option('--solver', type=click.Choice(['cb_orig', 'cb_gpu', 'cms', 'libfes', 'magma', 'mq', 'wdsat', 'xl'], case_sensitive=False), help='the external solver to be used')
@ click.option('--gen_only', default=False, is_flag=True, help='only generate equation systems')
@ click.option('--solve_only', default=False, is_flag=True, help='only solve an existing system and check solutions')
@ click.option('--log_path', default=Path("log.txt"), help='path to a log with the solution')
@ click.option('--inner_hybridation', '-h', default="-1", help='the number of variable that are not guessed in MQ', type=int)
@ click.option('--verbose', '-v', default=False, is_flag=True, help='control the output verbosity')
@ click.option('--seed', '-s', default=0, help='the seed for randomness replication', type=int)
@ click.option('--precompiled', default=False, is_flag=True, help='indicates if all relevant solvers are already compiled w.r.t. the parameters')
def main(q, n, m, o2, solver, gen_only, solve_only, log_path, inner_hybridation, verbose, seed, precompiled):
    if m == 0:
        m = 2*o2
    if n == 0:
        n = 3*o2
    set_random_seed(seed)
    if verbose:
        print(f"Generating Rainbow instance for seed={seed}, q={q}, m={m}, n={n}, o2={o2}...")
    rainbow = Rainbow(q, m, n, o2, seed=seed)
    system_folder_path, base_system_name = declare_paths(seed, q, o2, m, n)
    setup_path = Path(system_folder_path, base_system_name + '.stp')
    solution_path = Path(system_folder_path, base_system_name + '.sol')
    M, N = compute_system_size(q, m, n)

    if not solve_only:
        # get the attack equations
        equations, solution = rainbow.differential_attack(verbose=verbose)
        EqSys = EquationSystem(equations, seed=seed, verbose=verbose, solution=solution[1:-1])
        assert (EqSys.M, EqSys.N) == (M, N)
        # save everything
        EqSys.save_all(system_folder_path, base_system_name)
        save_setup(rainbow, setup_path, verbose=verbose)
    else:
        if verbose:
            print("Skipping the attack equations generation...")

    if gen_only or solver == None:
        exit()

    equations_path = get_eq_path(seed, q, o2, m, n, M, N, solver)
    invoke_solver(solver, equations_path, q, M, N, log_path=log_path,
                  inner_hybridation=inner_hybridation, precompiled=precompiled)

    if not solve_only:
        sol_path = get_sol_path(seed, q, o2, m, n, M, N, weil=use_weil(solver))
        check_cmd = f"sage check_solution.sage --q {q} --sol_path {sol_path} --solver {solver}"
        proc = sp.run(check_cmd, stdout=sp.PIPE, stderr=sp.STDOUT, shell=True)
        print(proc.stdout.decode())


if __name__ == '__main__':
    main()
