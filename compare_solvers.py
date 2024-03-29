#!/usr/bin/env python3

from functools import reduce
from operator import itemgetter
from pathlib import Path
from prettytable import PrettyTable
import click
import datetime
import json
import statistics
import subprocess as sp
from compile_solver import compile_solver
from invoke_solver import invoke_solver
from config_utils import defaults, get_eq_path, get_sol_path, solvers_to_skip, use_weil, get_rainbow_dims, get_weil_dims


def sec_to_str(t):
    units = reduce(lambda ll, b: divmod(ll[0], b) + ll[1:], [(t*1000,), 1000, 60, 60])
    return "{:01.0f}:{:02.0f}:{:02.0f}.{:02.0f}".format(*units[:-1], *[x/10 for x in units[-1:]])


@ click.command()
@ click.option('--q', default=2, type=click.Choice(['2', '16']), help='field size')
@ click.option('--o2_min', default=10, type=int, help='lower bound for o2', )
@ click.option('--o2_max', default=16, type=int, help='upper bound for o2')
@ click.option('--iterations', default=2, type=int, help='number of iterations for each parameter set')
@ click.option('--log_path_brief', default=defaults("comparison_brief"), help='the path to the brief log')
@ click.option('--log_path_verbose', default=defaults("comparison_verbose"), help='the path to the verbose log')
@ click.option('--to_skip', '-s', default=solvers_to_skip(), type=click.Choice(defaults("solvers"), case_sensitive=False), multiple=True, help='the solvers to be skipped')
@ click.option('--timeout', '-t', default=defaults("timeout"),  help='the maximum time (in seconds) allowed for running the solver')
def main(q, o2_min, o2_max, iterations, log_path_brief, log_path_verbose, to_skip, timeout):

    # Set up main parameters
    solvers = [s for s in defaults("solvers") if s not in set(to_skip)]
    q = int(q)
    o2_range = range(o2_min, o2_max + 1, 2)
    outcomes = ["success", "failure"]

    # Set up logging
    json_path = Path("comparison.json")
    with open(json_path, "w") as j:
        results = dict.fromkeys([f"q={str(q)}"], dict.fromkeys(["o2=" + str(o2)
                                for o2 in o2_range], dict.fromkeys(solvers, None)))
        json.dump(results, j)

    def print_and_log(to_log, to_print=None, include_brief=False):
        if to_print == None:
            to_print = to_log
        if not to_print == "":
            print(to_print)
        # Write only into the verbose log by default
        log_paths = [log_path_brief, log_path_verbose][not include_brief:]
        for path in log_paths:
            with open(path, 'a') as f:
                f.write(to_log + "\n")

    # Set up formatting
    star_length = 105
    stars = '*' * star_length

    # Set up tables
    T = PrettyTable()
    res_col_name = f"Successes (out of {iterations})"
    T.field_names = ["Solver", f"{res_col_name}",
                     "Avg time (s)", "Stdev of time (s)", "Avg resident memory (MB)"]
    T.align = "c"

    # Use colors (from https://stackoverflow.com/questions/43583847/python-pretty-table-with-color-output)
    red = "\033[0;31;40m"
    green = "\033[0;32;40m"
    reset = "\033[0m"
    colors = [green, red]

    # Sort the tables by results and times
    def sort_key(r):
        res, time = itemgetter(2, 3)(r)
        if reset in res:
            res = res.split(reset)[0]
            for color in colors:
                if color in res:
                    res = res.split(color)[1]
        return -int(res), time

    T.sortby = res_col_name
    T.sort_key = sort_key
    T_color = T.copy()

    # Start the comparison
    print_and_log(
        f"{stars}\nStarting solver comparison, check {log_path_verbose} for the detailed results.", include_brief=True)
    print_and_log(f"Current datetime: {datetime.datetime.now().isoformat(' ', 'seconds')}", include_brief=True)
    for o2 in o2_range:
        m, n, M, N = get_rainbow_dims(o2)
        gen_msg = f"Generating equations for q = {q: 2}, o2 = {o2: 2}  ==>  M = {M: 2}, N = {N: 2} (before Weil descent); {iterations: 2} iterations..."
        print_and_log(f"\n\n{stars}\n{gen_msg}\n{stars}", include_brief=True)
        solver_stats = {solver: {"successes": 0, "times": [], "memories": [],
                                 "mean_time": None, "stdev_time": None, "mean_memory": None} for solver in solvers}

        # Go through all iterations for the given parameters
        for seed in range(iterations):
            try:
                gen_cmd = f"sage rainbow_attacks.sage --seed {seed} --q {q} --o2 {o2} --m {m} --n {n} --gen_only"
            except Exception as e:
                print("An error ocurred during generating a Rainbow instance:", e)
                continue

            sp.call(gen_cmd, shell=True)
            for solver in solvers:
                # Determine if Weil descent should be used
                weil = use_weil(solver)
                if weil:
                    M_weil, N_weil = get_weil_dims(q, M, N)
                _q, _M, _N = (2, M_weil, N_weil) if weil else (q, M, N)

                # Compile the solver for each parameter set if needed
                if seed == 0 and solver in defaults("solvers_to_compile"):
                    try:
                        out = compile_solver(solver, _q, _M, _N)
                        print_and_log(out, to_print="")
                    except Exception as e:
                        print_and_log(str(e), to_print="")

                print_and_log(
                    f"\n{stars}\nCurrent datetime: {datetime.datetime.now().isoformat(' ', 'seconds')}", to_print="")
                print_and_log(f"Solver: {solver}\n", to_print="")

                # Measure the time and memory usage of the active process and all its subprocesses
                try:
                    eq_path = get_eq_path(seed, q, o2, m, n, M, N, solver)
                    out, time_taken, rss = invoke_solver(solver, eq_path, _q, _M, _N, precompiled=True, timeout=timeout)
                    print_and_log(out, to_print="")
                    sol_path = get_sol_path(seed, q, o2, m, n, M, N, weil)
                    check_cmd = f"sage check_solution.sage --q {q} --sol_path {sol_path} --solver {solver}"
                    proc = sp.run(check_cmd, stdout=sp.PIPE, stderr=sp.STDOUT, shell=True)
                    print_and_log(proc.stdout.decode(), to_print="")
                    # successful solutions should yield code 0
                    code = int('Attack successful!' not in proc.stdout.decode())
                except Exception as e:
                    print_and_log(str(e), to_print="")
                    time_taken = 0
                    rss = 0
                    code = 1

                # Save the iteration results
                solver_stats[solver]["successes"] += (1-code)
                solver_stats[solver]["times"].append(time_taken)
                solver_stats[solver]["memories"].append(rss)
                print_and_log(f"Result: {outcomes[code]}", to_print="")
                print_and_log(f"Time:   {sec_to_str(time_taken)}", to_print="")
                print_and_log(f"Memory: {( rss / 1000000): .2f} MB", to_print="")
                print_and_log(stars, to_print="")

        # Save the aggregated results
        for solver in solvers:
            successes = str(solver_stats[solver]["successes"])
            solver_stats[solver]["successes"] = f"{successes} of {iterations}"
            mean_time = statistics.mean(solver_stats[solver]["times"])
            solver_stats[solver]["mean_time"] = round(mean_time, 2)
            stdev_time = round(statistics.stdev(solver_stats[solver]["times"]), 2)
            solver_stats[solver]["stdev_time"] = stdev_time
            mean_memory = statistics.mean(solver_stats[solver]["memories"])
            mean_memory = round(mean_memory / 1000000, 2)
            solver_stats[solver]["mean_memory"] = mean_memory
            del solver_stats[solver]['times']
            del solver_stats[solver]['memories']

            T.add_row([solver, successes, sec_to_str(mean_time), stdev_time, f"{mean_memory:.1f}"])
            colored_successes = colors[successes != str(iterations)] + successes + reset
            T_color.add_row([solver, colored_successes, sec_to_str(mean_time), stdev_time, f"{mean_memory: .1f}"])

            # Update the JSON with results
            with open(json_path) as j:
                results = json.load(j)
                results[f"q={q}"][f"o2={o2}"][solver] = solver_stats[solver]
            with open(json_path, 'w') as j:
                json.dump(results, j)

        # Show the result table
        print_and_log(T.get_string(), to_print=T_color.get_string(), include_brief=True)
        T.clear_rows()
        T_color.clear_rows()


if __name__ == '__main__':
    main()
