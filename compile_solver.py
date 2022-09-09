from pathlib import Path
from subprocess import Popen
import click
import json


def create_wdsat_config(wdsat_path, M, N):
    with open(Path(wdsat_path, "src", "config.h"), 'w') as file:
        file.write("""
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
# define __MAX_XEQ_SIZE__ {4}""".format(N, M, N + 1, int(N*(N-1)/2), N * (N + 1), int(M*(M-1)/2)))


def check_params(status_path, q, M, N):
    if status_path.exists():
        with open(status_path, 'r') as f:
            params = json.load(f)
            if params['q'] == q and params['M'] == M and params['N'] == N:
                return True
    return False


@ click.command()
@ click.option('--solver', type=click.Choice(['xl', 'wdsat'], case_sensitive=False), help='the external solver to be compiled')
@ click.option('--q', help='field characteristic - needed for XL compilation', type=int)
@ click.option('--m', help='number of equations - needed for XL and WDSAT compilation', type=int)
@ click.option('--n', help='number of variables - needed for XL and WDSAT compilation', type=int)
@ click.option('--xl_path', default=Path("..", "xl"), help='the path the XL solver folder: http://polycephaly.org/projects/xl', type=str)
@ click.option('--wdsat_path', default=Path("..", "WDSat"), help='the path the WDSat solver folder: https://github.com/mtrimoska/WDSat', type=str)
def main(solver, q, m, n, xl_path, wdsat_path,):
    if not solver:
        print("Please specify a solver.")
        exit()

    M = m
    N = n

    if solver == 'xl':
        xl_status_path = Path("xl_status.json")
        compiled = check_params(xl_status_path, q, M, N)
        if compiled and Path(xl_path, "xl").exists():
            print("\nThe XL solver is already compiled.")
        else:
            make_cmd = f"make -C {str(xl_path)} clean && make -C {str(xl_path)} Q={q} M={M} N={N} -Wno-unused-result -Wno-class-memaccess"
            print("\nCompiling the XL solver...")
            Popen(make_cmd, shell=True).wait()
            with open(xl_status_path, 'w') as f:
                params = {'q': q, 'M': M, 'N': N}
                json.dump(params, f)

    if solver == 'wdsat':
        wdsat_status_path = Path("wdsat_status.json")
        compiled = check_params(wdsat_status_path, q, m, n)
        if compiled and Path(wdsat_path, "wdsat_solver").exists():
            print("\nThe WDSat solver is already compiled.")
        else:
            create_wdsat_config(wdsat_path, m, n)
            suppressor_flag = ""
            src_path = str(Path(wdsat_path, "src"))
            make_cmd = f"make -C {src_path} {suppressor_flag}"
            print("\nCompiling the WDSat solver...")
            Popen(make_cmd, shell=True).wait()
            with open(wdsat_status_path, 'w') as f:
                params = {'q': q, 'M': M, 'N': N}
                json.dump(params, f)


if __name__ == '__main__':
    main()
