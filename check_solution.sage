#!/usr/bin/env sage
import click
import re
from config_utils import defaults, get_log_format
load("equation_utils.sage")


def load_solution(solution_path, q, ext_format='alg'):
    with open(solution_path, 'r') as f:
        if ZZ(q).is_prime():
            sol_str = f.readlines()[0].strip()
            return list(sage_eval(sol_str))
        elif ext_format == 'alg':
            sol_strs = f.readlines()[0].strip().split("(")[1].split(")")[0].split(",")
            return list(GF(q)(sol_str) for sol_str in sol_strs)
        elif ext_format == 'hex':
            sol_strs = f.readlines()[0].strip().split("[")[1].split("]")[0].split(",")
            return list(str_to_elt(q, sol_str) for sol_str in sol_strs)
        elif ext_format == 'hex_to_weil':
            # used when the solution file is hex encoded, but we need a Weil format (e.g., for Fukuoka challenge answers)
            sol_strs = f.readlines()[0].strip().split("[")[1].split("]")[0].split(",")
            hex_sols = list(str_to_elt(q, sol_str) for sol_str in sol_strs)
            return([bit for hex_sol in hex_sols for bit in weil_decomposition(hex_sol)])
        else:
            raise Exception("ext_format must be either 'alg', 'hex', or 'hex_to_weil'")


def get_solution_from_log(log_path, format, N, q, ext_deg=1):
    assert format in ['cb', 'cms', 'magma', 'mq', 'wdsat', 'xl']
    with open(log_path, 'r') as f:
        z = GF(q).gens()[0]
        zs = [z ^ i for i in range(ext_deg)]
        variables = []
        found = False
        sol_str = ""
        for line in f.readlines():
            if format == 'cms':
                Nw = N * ext_deg
                if line[0] == 'v':
                    variables += line.strip().split(" ")[1:]
                    if len(variables) < Nw:
                        continue
                    sol = [int((sgn(int(b)) + 1) / 2) for b in variables[:Nw]]
                    parts = [sol[ext_deg * i:ext_deg * i + ext_deg]
                             for i in range(len(sol) / ext_deg)]
                    return [linear_combination(bits, zs) for bits in parts]

            if format == 'cb':
                if "solution found: " in line:
                    found = True
                    continue
                if found:
                    sol = sage_eval(line.strip())
                    parts = [sol[ext_deg * i:ext_deg * i + ext_deg]
                             for i in range(len(sol) / ext_deg)]
                    return [linear_combination(bits, zs) for bits in parts]

            if format == 'magma':
                # find solutions even across multiple lines and take the first one
                if line[:3] == "[ <":
                    found = True
                if found:
                    sol_str += line
                    if ">" in line:
                        end = sol_str.index(">")
                        sol = [GF(q)(x) for x in sol_str[3:end].split(",")]
                        parts = [sol[ext_deg * i:ext_deg * i + ext_deg]
                                 for i in range(len(sol) / ext_deg)]
                        return [linear_combination(bits, zs) for bits in parts]

            if format == 'mq':
                if "solution found : " in line:
                    sol = [int(b) for b in line.split(
                        "solution found : ")[1][1:-2].split(", ")]
                    parts = [sol[ext_deg * i:ext_deg * i + ext_deg]
                             for i in range(len(sol) / ext_deg)]
                    return [linear_combination(bits, zs) for bits in parts]

            if format == 'wdsat':
                line = line.strip()
                if line == "":
                    continue
                if re.match('^[0-1]*$', line):
                    sol = [ZZ(b) for b in list(line)]
                    parts = [sol[ext_deg * i:ext_deg * i + ext_deg]
                             for i in range(len(sol) / ext_deg)]
                    return [linear_combination(bits, zs) for bits in parts]

            if format == 'xl':
                if "  is sol" in line:
                    sol = line.split("  is sol")[0].split(" ")
                    return [str_to_elt(q, c) for c in sol]
    return None


@ click.command()
@ click.option('--q', type=int, help='the size of the finite field')
@ click.option('--sol_path', help='path to a log with the expected solution')
@ click.option('--log_path', default=defaults("log_path"), help='path to a log with the found solution')
@ click.option('--solver', type=click.Choice(defaults("solvers"), case_sensitive=False), help='the solver used to find the solution')
@ click.option('--ext_format', default='alg', type=click.Choice(['alg', 'hex', 'hex_to_weil'], case_sensitive=False), help='the format of the expected solution')
def main(q, sol_path, log_path, solver, ext_format):
    # load the expected solution
    try:
        if q is None:
            raise Exception("No field order q specified.")
        if sol_path is None:
            raise Exception("No solution path specified.")
        solution_expected = load_solution(sol_path, q, ext_format)
    except Exception as e:
        print("An error ocurred during loading the solution: ", e)
        return False
    N = len(solution_expected)

    # load the actual solution from a log
    try:
        if solver is None:
            raise Exception("No solver specified.")
        log_format = get_log_format(solver)
        solution_found = list(get_solution_from_log(log_path, format=log_format, N=N, q=q))
    except Exception as e:
        print("An error ocurred during parsing the log: ", e)
        solution_found = None

    # handle potential dummy variables
    try:
        if solver in ["libfes", "mq"]:
            solution_found = solution_found[:N]
    except TypeError:
        pass

    # check if the two solutions agree
    print(f"\n{'First solution found: ' : <25} {solution_found} ")
    print(f"{'Expected solution: ' : <25} {solution_expected}\n")
    success = solution_found == solution_expected
    if success:
        print("Attack successful!\n")
    else:
        print("Attack NOT successful. :(\n")
    return success


if __name__ == '__main__':
    main()
