#!/usr/bin/env python3

from pathlib import Path
import subprocess as sp
import click
import json
from config_utils import defaults


def compile_solver(solver, q, m, n, cb_orig_path=defaults("cb_orig_path"), wdsat_path=defaults("wdsat_path"), xl_path=defaults("xl_path")):
    if not solver:
        raise Exception("No solver specified.")

    M = m
    N = n

    Path("tmp").mkdir(exist_ok=True)

    if solver == 'cb_orig':
        cb_orig_status_path = Path("tmp", "cb_orig_status.json")
        compiled = check_params(cb_orig_status_path, q, m, n)
        if compiled and Path(cb_orig_path, "LinBlockLanczos").exists() and Path(cb_orig_path, "CheckCandidate").exists():
            out = "The crossbred (original) solver is already compiled."
        else:
            suppressor_flag = "-Wno-unused-result -Wno-format -Wno-shift-count-overflow"
            src = Path(cb_orig_path, "LinBlockLanczos.c")
            binary = Path(cb_orig_path, "LinBlockLanczos")
            if binary.exists():
                binary.unlink()
            trunc_var = get_trunc_var(M, N)
            params = f"-DNBVARS={N} -DTRUNC_VAR={trunc_var} -DPRINTVARS={N} -DNBPOLS={M}"
            gcc_cmd = f"gcc {params} -Ofast -march=native {suppressor_flag} -o {binary} {src}"
            print("Compiling the crossbred (original) solver (linear algebra)...")
            out = sp.run(gcc_cmd, stdout=sp.PIPE, stderr=sp.STDOUT, shell=True).stdout.decode()

            src = Path(cb_orig_path, "CheckCandidates.c")
            binary = Path(cb_orig_path, "CheckCandidates")
            if binary.exists():
                binary.unlink()
            gcc_cmd = f"gcc {params} -o {binary} {src}"
            print("Compiling the crossbred (original) solver (checking)...")
            out += sp.run(gcc_cmd, stdout=sp.PIPE, stderr=sp.STDOUT, shell=True).stdout.decode()

            with open(cb_orig_status_path, 'w') as f:
                params = {'q': int(q), 'M': int(M), 'N': int(N)}
                json.dump(params, f)

    if solver == 'wdsat':
        wdsat_status_path = Path("tmp", "wdsat_status.json")
        compiled = check_params(wdsat_status_path, q, m, n)
        if compiled and Path(wdsat_path, "wdsat_solver").exists():
            out = "The WDSat solver is already compiled."
        else:
            create_wdsat_config(wdsat_path, m, n)
            suppressor_flag = ""
            src_path = str(Path(wdsat_path, "src"))
            make_cmd = f"make -C {src_path} clean && make -C {src_path} {suppressor_flag}"
            print("Compiling the WDSat solver...")
            out = sp.run(make_cmd, stdout=sp.PIPE, stderr=sp.STDOUT, shell=True).stdout.decode()
        with open(wdsat_status_path, 'w') as f:
            params = {'q': int(q), 'M': int(M), 'N': int(N)}
            json.dump(params, f)

    if solver == 'xl':
        xl_status_path = Path("tmp", "xl_status.json")
        compiled = check_params(xl_status_path, q, M, N)
        if compiled and Path(xl_path, "xl").exists():
            out = "The XL solver is already compiled."
        else:
            make_cmd = f"make -C {str(xl_path)} clean && make -C {str(xl_path)} Q={q} M={M} N={N} -Wno-unused-result -Wno-class-memaccess"
            print("Compiling the XL solver...")
            out = sp.run(make_cmd, stdout=sp.PIPE, stderr=sp.STDOUT, shell=True).stdout.decode()
            with open(xl_status_path, 'w') as f:
                params = {'q': int(q), 'M': int(M), 'N': int(N)}
                json.dump(params, f)
    return out


def create_wdsat_config(wdsat_path, M, N):
    with open(Path(wdsat_path, "src", "config.h"), 'w') as f:
        f.write("""
//enable the XG-ext module (must use anf as input)
# define __XG_ENHANCED__

//find all solutions instead of only one
# define __FIND_ALL_SOLUTIONS__

/** Rainbow : N={0} M={1} **/
# ifdef __XG_ENHANCED__
# define __MAX_ANF_ID__ {2} // make it +1
# define __MAX_DEGREE__ 3 // make it +1
# endif
# define __MAX_ID__ {3}
# define __MAX_BUFFER_SIZE__ 200000
# define __MAX_EQ__ {5}
# define __MAX_EQ_SIZE__ 4 //make it +1
# define __MAX_XEQ__ {1}
# define __MAX_XEQ_SIZE__ {4}""".format(N, M, N + 1, int(N*(N+1)/2), N * (N + 1), int(M*(M-1)/2)))


def check_params(status_path, q, M, N):
    if status_path.exists():
        with open(status_path, 'r') as f:
            try:
                params = json.load(f)
                if params['q'] == q and params['M'] == M and params['N'] == N:
                    return True
            except (KeyError, json.decoder.JSONDecodeError):
                pass
    return False


def get_trunc_var(M, N):
    # preliminary naive version
    return round(3 * N / 4)


@ click.command()
@ click.option('--solver', type=click.Choice(defaults("solvers_to_compile"), case_sensitive=False), help='the external solver to be compiled')
@ click.option('--q', help='field characteristic - needed for XL compilation', type=int)
@ click.option('--m', help='number of equations - needed for XL and WDSAT compilation', type=int)
@ click.option('--n', help='number of variables - needed for XL and WDSAT compilation', type=int)
@ click.option('--cb_orig_path', default=defaults("cb_orig_path"), help='the path the crossbred (original) solver folder', type=str)
@ click.option('--wdsat_path', default=defaults("wdsat_path"), help='the path the WDSat solver folder: https://github.com/mtrimoska/WDSat', type=str)
@ click.option('--xl_path', default=defaults("xl_path"), help='the path the XL solver folder: http://polycephaly.org/projects/xl', type=str)
def main(solver, q, m, n, cb_orig_path, wdsat_path, xl_path):
    try:
        compile_solver(solver, q, m, n, cb_orig_path, wdsat_path, xl_path)
    except Exception as e:
        print("An error ocurred during compiling a solver: ", e)


if __name__ == '__main__':
    main()
