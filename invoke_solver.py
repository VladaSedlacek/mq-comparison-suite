#!/usr/bin/env python3

from pathlib import Path
import click
import subprocess as sp
import psutil
import time
from compile_solver import compile_solver


def invoke_solver(solver, equations_path, q, m, n, log_path=Path(".", "log.txt"), cb_gpu_path=Path("..", "mqsolver"), cb_orig_path=Path("..", "crossbred"), cms_path=Path("..", "cryptominisat", "build"), libfes_path=Path("..", "libfes-lite", "build"), magma_path=Path("magma"), mq_path=Path("..", "mq"), wdsat_path=Path("..", "WDSat"), xl_path=Path("..", "xl"), inner_hybridation=-1, precompiled=False, timeout=1000):

    if not solver:
        print("Please specify a solver.")
        exit()

    if not Path(equations_path).exists():
        print("Please specify an existing equation file.")
        exit()

    if solver == 'cb_orig':
        linalg_path = Path(cb_orig_path, "LinBlockLanczos")
        check_path = Path(cb_orig_path, "CheckCandidates")
        if not (precompiled and linalg_path.exists() and check_path.exists()):
            compile_solver('cb_orig', q, m, n, cb_orig_path)
        print("Starting the crossbred (original) solver...")
        solve_cmd = f"{linalg_path} {equations_path} | tee {log_path}"
        start_time = time.time()
        proc = sp.Popen(solve_cmd, stdout=sp.PIPE, stderr=sp.STDOUT, shell=True)
        while proc.poll() is None:
            rss = psutil.Process(proc.pid).memory_info().rss
            # the memory used in the checking part is neglected
            proc.wait(timeout)
        out = proc.communicate()[0].decode()
        candidates = [cand for cand in out.strip().split("\n") if "@" in cand]
        out += "\n"
        for cand in candidates:
            out += f"Candidate: {cand}\n"
            check_cmd = f"echo {cand} | {check_path} {equations_path}"
            check_out = sp.run(check_cmd, stdout=sp.PIPE, stderr=sp.STDOUT, shell=True).stdout.decode()
            res = check_out.split("\n")
            # ensure compatibility with mqsolver log
            if "solution found :)" in res:
                assert res[1] == '0' * m
                out += f"solution found: \n[{res[0]}]\n"
            else:
                out += f"does not work: \n[{res[0]}]\n\tevaluates to {res[1]}"
        time_taken = time.time() - start_time

    else:
        cwd = Path.cwd()

        if solver == 'cms':
            print("Starting the CryptoMiniSat solver...")
            p = Path(cms_path, "cryptominisat5")
            solve_cmd = f"{p} --verb 0 {equations_path}"

        if solver == 'cb_gpu':
            print("Starting the crossbred (GPU) solver...")
            solve_cmd = f"./solve.py -d 3 -k 16 -t 20 -v {Path(cwd, equations_path)}"
            cwd = Path(cb_gpu_path)

        if solver == 'libfes':
            print("Starting the libfes solver...")
            p = Path(libfes_path, "benchmark", "demo")
            solve_cmd = f"{p} < {equations_path}"

        if solver == 'magma':
            print("Starting Magma...")
            solve_cmd = f"{magma_path} < {equations_path}"

        if solver == 'mq':
            print("Starting the MQ solver...")
            inner_hybridation_arg = " --inner-hybridation " + str(inner_hybridation) if inner_hybridation != -1 else ""
            # Use the non-vectorized version for less than 8 variables
            if 3 <= n and n <= 8:
                # the first ineqality seems to prevent an infinite loop
                binary = "monica"
            else:
                binary = "monica_vector"
            p = Path(mq_path, f"{binary}")
            solve_cmd = f"{p}{inner_hybridation_arg} < {equations_path}"

        if solver == 'wdsat':
            if not (precompiled and Path(wdsat_path, "wdsat_solver").exists()):
                compile_solver('wdsat', q, m, n, wdsat_path)
            print("Starting the WDSat solver...")
            p = Path(wdsat_path, "wdsat_solver")
            solve_cmd = f"{p} -i {equations_path}"

        if solver == 'xl':
            if not (precompiled and Path(xl_path, "xl").exists()):
                compile_solver('xl', q, m, n, xl_path)
            print("Starting the XL solver...")
            p = Path(xl_path, "xl")
            solve_cmd = f"{p} --challenge {equations_path} --all"

        start_time = time.time()
        proc = sp.Popen(solve_cmd, stdout=sp.PIPE, stderr=sp.STDOUT, shell=True, cwd=cwd)
        while proc.poll() is None:
            rss = psutil.Process(proc.pid).memory_info().rss
            proc.wait(timeout)
        time_taken = time.time() - start_time
        out = proc.communicate()[0].decode()

    with open(log_path, 'w') as f:
        f.write(out)

    return out, time_taken, rss


@ click.command()
@ click.option('--solver', type=click.Choice(['cb_gpu', 'cb_orig', 'cms', 'libfes', 'magma', 'mq', 'wdsat', 'xl'], case_sensitive=False), help='the external solver to be used')
@ click.option('--equations_path', '-e', help='the path to the equation system', type=str)
@ click.option('--q', help='field characteristic - needed for XL compilation', type=int)
@ click.option('--m', help='number of equations - needed for XL and WDSAT compilation', type=int)
@ click.option('--n', help='number of variables - needed for XL and WDSAT compilation', type=int)
@ click.option('--log_path', '-l', default=Path(".", "log.txt"), help='the path to the output log', type=str)
@ click.option('--cb_gpu_path', default=Path("..", "mqsolver"), help='the path the crossbred solver folder: https://github.com/kcning/mqsolver', type=str)
@ click.option('--cb_orig_path', default=Path("..", "crossbred"), help='the path the crossbred (original) solver folder', type=str)
@ click.option('--cms_path', default=Path("..", "cryptominisat", "build"), help='the path the CMS solver folder: https://github.com/msoos/cryptominisat', type=str)
@ click.option('--libfes_path', default=Path("..", "libfes-lite", "build"), help='the path the libfes solver folder: https://github.com/cbouilla/libfes-lite', type=str)
@ click.option('--magma_path', default=Path("magma"), help='the path the Magma binary: https://magma.maths.usyd.edu.au', type=str)
@ click.option('--mq_path', default=Path("..", "mq"), help='the path the MQ solver folder: https://gitlab.lip6.fr/almasty/mq', type=str)
@ click.option('--wdsat_path', default=Path("..", "WDSat"), help='the path the WDSat solver folder: https://github.com/mtrimoska/WDSat', type=str)
@ click.option('--xl_path', default=Path("..", "xl"), help='the path the XL solver folder: http://polycephaly.org/projects/xl', type=str)
@ click.option('--inner_hybridation', '-h', default="-1", help='the number of variable that are not guessed in MQ', type=int)
@ click.option('--precompiled', default=False, is_flag=True, help='indicates if all relevant solvers are already compiled w.r.t. the parameters')
@ click.option('--timeout', default=1000,  help='the maximum time allowed for running the solver')
def main(solver, equations_path, q, m, n, log_path, cb_gpu_path, cb_orig_path, cms_path, libfes_path, magma_path, mq_path, wdsat_path, xl_path, inner_hybridation, precompiled, timeout):
    out, time_taken, rss = invoke_solver(solver, equations_path, q, m, n, log_path, cb_gpu_path, cb_orig_path,
                                         cms_path, libfes_path, magma_path, mq_path, wdsat_path, xl_path, inner_hybridation, precompiled, timeout)
    print(out)
    print(f"Time taken: {time_taken: .2f} s")
    print(f"Resident memory used: {( rss / 1000000): .2f} MB")


if __name__ == '__main__':
    main()
