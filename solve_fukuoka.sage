#!/usr/bin/env sage

import click
from pathlib import Path
import subprocess as sp
from config_utils import defaults, get_eq_format, use_weil
from invoke_solver import invoke_solver

load("equation_utils.sage")


@ click.command()
@ click.option('--fukuoka_challenge_path', help='path to the equations file downloaded from https://www.mqchallenge.org/')
@ click.option('--fukuoka_answer_path', help='path to the solution file')
@ click.option('--solver', type=click.Choice(defaults("solvers"), case_sensitive=False), help='the external solver to be used')
@ click.option('--log_path', default=defaults("log_path"), help='path to a log with the solution')
def main(fukuoka_challenge_path, fukuoka_answer_path, solver, log_path):
    if fukuoka_answer_path == None:
        fukuoka_answer_path = fukuoka_challenge_path + "-answer"

    try:
        folder = Path(fukuoka_challenge_path).parent
        base_path = Path(fukuoka_challenge_path).stem
        EqSys = load_fukuoka(fukuoka_challenge_path)
        EqSys.save_all(folder, base_path, save_sol=False)
    except Exception as e:
        print("An error ocurred during converting a Fukuoka challenge: ", e)
    if solver == None:
        exit()

    # apply Weil descent if needed
    if use_weil(solver) and EqSys.weil is not None:
        q, M, N = (EqSys.weil.q, EqSys.weil.M, EqSys.weil.N)
    else:
        q, M, N = (EqSys.q, EqSys.M, EqSys.N)

    try:
        equations_path = Path(folder, f"{base_path}_M_{M}_N_{N}.{get_eq_format(solver)}")
        invoke_solver(solver, equations_path, q, M, N, log_path=log_path)
    except Exception as e:
        print("An error ocurred during invoking a solver: ", e)
    check_cmd = f"sage check_solution.sage --q {q} --sol_path {fukuoka_answer_path} --solver {solver}"
    proc = sp.run(check_cmd, stdout=sp.PIPE, stderr=sp.STDOUT, shell=True)
    print(proc.stdout.decode())


if __name__ == '__main__':
    main()
