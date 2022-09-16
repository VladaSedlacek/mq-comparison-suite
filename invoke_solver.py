#!/usr/bin/env python3

from pathlib import Path
from subprocess import Popen
import click
import os


@ click.command()
@ click.option('--solver', type=click.Choice(['cms', 'crossbred', 'libfes', 'magma', 'mq', 'wdsat', 'xl'], case_sensitive=False), help='the external solver to be used')
@ click.option('--equations_path', '-e', help='the path to the equation system', type=str)
@ click.option('--q', help='field characteristic - needed for XL compilation', type=int)
@ click.option('--m', help='number of equations - needed for XL and WDSAT compilation', type=int)
@ click.option('--n', help='number of variables - needed for XL and WDSAT compilation', type=int)
@ click.option('--log_path', '-l', default=Path(".", "log.txt"), help='the path to the output log', type=str)
@ click.option('--cms_path', default=Path("..", "cryptominisat", "build"), help='the path the CMS solver folder: https://github.com/msoos/cryptominisat', type=str)
@ click.option('--crossbred_path', default=Path("..", "mqsolver"), help='the path the crossbred solver folder: https://github.com/kcning/mqsolver', type=str)
@ click.option('--libfes_path', default=Path("..", "libfes-lite", "build"), help='the path the libfes solver folder: https://github.com/cbouilla/libfes-lite', type=str)
@ click.option('--magma_path', default=Path("magma"), help='the path the Magma binary: https://magma.maths.usyd.edu.au', type=str)
@ click.option('--mq_path', default=Path("..", "mq"), help='the path the MQ solver folder: https://gitlab.lip6.fr/almasty/mq', type=str)
@ click.option('--wdsat_path', default=Path("..", "WDSat"), help='the path the WDSat solver folder: https://github.com/mtrimoska/WDSat', type=str)
@ click.option('--xl_path', default=Path("..", "xl"), help='the path the XL solver folder: http://polycephaly.org/projects/xl', type=str)
@ click.option('--inner_hybridation', '-h', default="-1", help='the number of variable that are not guessed in MQ', type=int)
@ click.option('--precompiled', default=False, is_flag=True, help='indicates if all relevant solvers are already compiled w.r.t. the parameters')
def main(solver, equations_path, q, m, n, log_path, cms_path, crossbred_path, libfes_path, magma_path, mq_path, wdsat_path, xl_path, inner_hybridation, precompiled):
    if not solver:
        print("Please specify a solver.")
        exit()

    if not Path(equations_path).exists():
        print("Please specify an existing equation file.")
        exit()

    current_path = Path().cwd()

    if solver == 'cms':
        Popen(" > {}".format(str(log_path)), shell=True).wait()
        print("\nStarting the CryptoMiniSat solver...")
        cms_solve_cmd = "{} --verb 0 {}".format(
            str(Path(cms_path, "cryptominisat5")), str(equations_path) +
            " | tee " + str(log_path))
        Popen(cms_solve_cmd, shell=True).wait()

    if solver == 'crossbred':
        print("\nStarting the crossbred solver...")
        os.chdir(crossbred_path)
        crossbred_solve_cmd = "./solve.py -d 3 -k 16 -t 20 -v -o {} {}".format(
            str(Path(current_path, log_path)), str(Path(current_path, equations_path)))
        Popen(crossbred_solve_cmd, shell=True).wait()
        os.chdir(current_path)

    if solver == 'libfes':
        print("\nStarting the libfes solver...")
        mq_solve_cmd = "{} < {}".format(
            str(Path(libfes_path, "benchmark", "demo")), str(equations_path)) + " | tee " + str(log_path)
        Popen(mq_solve_cmd, shell=True).wait()

    if solver == 'magma':
        Popen(" > {}".format(str(log_path)), shell=True).wait()
        print("\nStarting Magma...")
        magma_solve_cmd = f"{magma_path} < {equations_path} | tee {log_path}"
        Popen(magma_solve_cmd, shell=True).wait()

    if solver == 'mq':
        print("\nStarting the MQ solver...")
        inner_hybridation_arg = " --inner-hybridation " + \
            str(inner_hybridation) if inner_hybridation != -1 else ""
        # Use the non-vectorized version for less than 8 variables
        if 3 <= n and n <= 8:
            # the first ineqality seems to prevent and infinite loop
            binary = "monica"
        else:
            binary = "monica_vector"
        mq_solve_cmd = "{}{} < {}".format(
            str(Path(mq_path, f"{binary}")), inner_hybridation_arg, str(equations_path)) + " | tee " + str(log_path)
        Popen(mq_solve_cmd, shell=True).wait()

    if solver == 'wdsat':
        if precompiled and Path(wdsat_path, "wdsat_solver").exists():
            print("\nThe WDSat solver is already compiled.")
        else:
            make_cmd = f"python3 compile_solver.py --solver wdsat --q {q} --m {m} --n {n} --wdsat_path {wdsat_path}"
            Popen(make_cmd, shell=True).wait()
        print("\nStarting the WDSat solver...")
        wdsat_solve_cmd = "{} -i {}".format(
            str(Path(wdsat_path, "wdsat_solver")), str(equations_path) +
            " | tee " + str(log_path))
        Popen(wdsat_solve_cmd, shell=True).wait()

    if solver == 'xl':
        if precompiled and Path(xl_path, "xl").exists():
            print("\nThe XL solver is already compiled.")
        else:
            make_cmd = f"python3 compile_solver.py --solver xl --q {q} --m {m} --n {n} --xl_path {xl_path}"
            Popen(make_cmd, shell=True).wait()
        print("\nStarting the XL solver...")
        xl_solve_cmd = "{} --challenge {} --all".format(
            str(Path(xl_path, "xl")), str(equations_path)) + " | tee " + str(log_path)
        Popen(xl_solve_cmd, shell=True).wait()


if __name__ == '__main__':
    main()
