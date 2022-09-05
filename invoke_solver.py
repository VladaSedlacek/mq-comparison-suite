from pathlib import Path
from subprocess import Popen
import click
import json
import os


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


@ click.command()
@ click.option('--solver', type=click.Choice(['xl', 'crossbred', 'mq', 'libfes', 'wdsat', 'cms'], case_sensitive=False), help='the external solver to be used')
@ click.option('-e', '--equations_path', help='the path to the equation system', type=str)
@ click.option('--q', help='field characteristic - needed for XL compilation', type=int)
@ click.option('--m', help='number of equations - needed for XL and WDSAT compilation', type=int)
@ click.option('--n', help='number of variables - needed for XL and WDSAT compilation', type=int)
@ click.option('-l', '--log_path', default=Path(".", "log.txt"), help='the path to the output log', type=str)
@ click.option('--xl_path', default=Path("..", "xl"), help='the path the XL solver folder: http://polycephaly.org/projects/xl', type=str)
@ click.option('--crossbred_path', default=Path("..", "mqsolver"), help='the path the crossbred solver folder: https://github.com/kcning/mqsolver/', type=str)
@ click.option('--mq_path', default=Path("..", "mq"), help='the path the MQ solver folder: https://gitlab.lip6.fr/almasty/mq', type=str)
@ click.option('--libfes_path', default=Path("..", "libfes-lite", "build"), help='the path the libfes solver folder: https://github.com/cbouilla/libfes-lite', type=str)
@ click.option('--wdsat_path', default=Path("..", "WDSat"), help='the path the WDSat solver folder: https://github.com/mtrimoska/WDSat', type=str)
@ click.option('--cms_path', default=Path("..", "cryptominisat", "build"), help='the path the CMS solver folder: https://github.com/mtrimoska/WDSat', type=str)
@ click.option('-h', '--inner_hybridation', default="-1", help='the number of variable that are not guessed in MQ', type=int)
def main(equations_path, log_path, q, m, n, xl_path, crossbred_path, mq_path, libfes_path, wdsat_path, cms_path, solver, inner_hybridation):
    if not solver:
        print("Please specify a solver.")
        exit()

    current_path = Path().cwd()

    def check_params(path, q, m, n):
        if path.exists():
            with open(path, 'r') as f:
                params = json.load(f)
                if params['q'] == q and params['m'] == m and params['n'] == n:
                    return True
        return False

    if solver == 'xl':
        xl_status_path = Path("xl_status.json")
        compiled = check_params(xl_status_path, q, m, n)
        if compiled:
            print("\nThe XL solver is already compiled.")
        else:
            make_command = "make -C {} Q={} M={} N={} -Wno-unused-result -Wno-class-memaccess".format(
                str(xl_path), str(q), str(m), str(n)) + " > " + str(log_path)
            print("\nCompiling the XL solver...")
            Popen(make_command, shell=True).wait()
            with open(xl_status_path, 'w') as f:
                params = {'q': q, 'm': m, 'n': n}
                json.dump(params, f)
        print("\nStarting the XL solver...")
        xl_solve_command = "{} --challenge {} --all".format(
            str(Path(xl_path, "xl")), str(equations_path)) + " | tee -a " + str(log_path)
        Popen(xl_solve_command, shell=True).wait()

    if solver == 'crossbred':
        print("\nStarting the crossbred solver...")
        os.chdir(crossbred_path)
        crossbred_solve_command = "./solve.py -d 3 -k 16 -t 20 -v -o {} {}".format(
            str(Path(current_path, log_path)), str(Path(current_path, equations_path)))
        Popen(crossbred_solve_command, shell=True).wait()
        os.chdir(current_path)

    if solver == 'mq':
        print("\nStarting the MQ solver...")
        inner_hybridation_arg = " --inner-hybridation " + \
            str(inner_hybridation) if inner_hybridation != -1 else ""
        mq_solve_command = "{}{} < {}".format(
            str(Path(mq_path, "monica_vector")), inner_hybridation_arg, str(equations_path)) + " | tee " + str(log_path)
        Popen(mq_solve_command, shell=True).wait()

    if solver == 'libfes':
        print("\nStarting the libfes solver...")
        mq_solve_command = "{} < {}".format(
            str(Path(libfes_path, "benchmark", "demo")), str(equations_path)) + " | tee " + str(log_path)
        Popen(mq_solve_command, shell=True).wait()

    if solver == 'wdsat':
        print("\nCompiling the WDSat solver...")
        create_wdsat_config(wdsat_path, m, n)
        suppressor_flag = "-n"
        make_command = "make -C {0} {1} > {2}".format(
            str(Path(wdsat_path, "src")), suppressor_flag, str(log_path))
        Popen(make_command, shell=True).wait()
        print("\nStarting the WDSat solver...")
        wdsat_solve_command = "{} -i {}".format(
            str(Path(wdsat_path, "wdsat_solver")), str(equations_path) +
            " | tee -a " + str(log_path))
        Popen(wdsat_solve_command, shell=True).wait()

    if solver == 'cms':
        Popen(" > {}".format(str(log_path)), shell=True).wait()
        print("\nStarting the CryptoMiniSat solver...")
        cms_solve_command = "{} --verb 0 {}".format(
            str(Path(cms_path, "cryptominisat5")), str(equations_path) +
            " | tee -a " + str(log_path))
        Popen(cms_solve_command, shell=True).wait()


if __name__ == '__main__':
    main()
