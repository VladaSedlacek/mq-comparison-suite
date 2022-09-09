from functools import reduce
from operator import itemgetter
from pathlib import Path
from subprocess import call, Popen
from prettytable import PrettyTable
import click
import datetime
import statistics
import time


def secondsToStr(t):
    return "%d:%02d:%02d.%03d" % \
        reduce(lambda ll, b: divmod(ll[0], b) + ll[1:],
               [(t*1000,), 1000, 60, 60])


@ click.command()
@ click.option('--o2_min', default=10, help='lower bound for o2', type=int)
@ click.option('--o2_max', default=16, help='upper bound for o2', type=int)
@ click.option('--iterations', default=2, help='number of iterations for each parameter set', type=int)
@ click.option('--log_path_brief', default=Path("comparison_log_brief.txt"), help='the path to the brief log')
@ click.option('--log_path_verbose', default=Path("comparison_log_verbose.txt"), help='the path to the verbose log')
def main(o2_min, o2_max, iterations, log_path_brief, log_path_verbose):

    # Set up main parameters
    solvers = ['xl', 'crossbred', 'mq', 'libfes', 'wdsat', 'cms']
    q_range = [2, 16]
    o2_range = range(o2_min, o2_max + 1, 2)
    results = ["success", "failure"]

    # Set up logging
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
    left_pad = ' ' * int((star_length - 72)/2)

    # Set up tables
    T = PrettyTable()
    res_col_name = f"Successes (out of {iterations})"
    T.field_names = ["Solver", f"{res_col_name}", "Average time", "Standard deviation of time"]
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
    for q in q_range:
        for o2 in o2_range:
            m = 2 * o2
            n = 3 * o2
            gen_msg = f"Generating equations for q = {q: 2}, o2 = {o2: 2}, m = {m: 2}, n = {n: 2}; {iterations: 2} iterations..."
            print_and_log(f"\n\n{stars}\n{left_pad}{gen_msg}\n{stars}", include_brief=True)
            solver_stats = {solver: {"successes": 0, "times": []} for solver in solvers}

            # Go through all iterations for the given parameters
            for seed in range(iterations):
                gen_cmd = f"sage rainbow_attacks.sage --seed {seed} --q {q} --o2 {o2} --m {m} --n {n}"
                call(gen_cmd, shell=True)
                for solver in solvers:

                    # Compile the solver for each parameter set if needed
                    if seed == 0 and solver in ["xl", "wdsat"]:
                        compile_cmd = f"python3 compile_solver.py --solver {solver} --q {q} --m {m-1} --n {n-m-2} >> {log_path_verbose} 2>&1"
                        call(compile_cmd, shell=True)

                    #  Measure the solving time and result
                    solve_cmd = f"sage rainbow_attacks.sage --seed {seed} --q {q} --o2 {o2} --m {m} --n {n} --solver {solver} --solve_only --precompiled"
                    print_and_log(f"\n{stars}\nExecuting: {solve_cmd}\n", to_print="")
                    solve_cmd += f" 2>> {log_path_verbose} | tee -a {str(log_path_verbose)}"
                    try:
                        start_time = time.time()
                        code = call(solve_cmd + " | grep 'Attack successful!' --quiet", shell=True)
                        time_taken = time.time() - start_time
                    except Exception as e:
                        print_and_log(str(e))
                        continue

                    # Save the iteration results
                    solver_stats[solver]["successes"] += (1-code)
                    solver_stats[solver]["times"].append(time_taken)
                    print_and_log(f"Solver: {solver}", to_print="")
                    print_and_log(f"Result: {results[code]}", to_print="")
                    print_and_log(f"Time:   {secondsToStr(time_taken)}", to_print="")
                    print_and_log(stars, to_print="")

            # Save the aggregated results
            for solver in solvers:
                successes = str(solver_stats[solver]["successes"])
                mean = secondsToStr(statistics.mean(solver_stats[solver]["times"]))
                stdev = secondsToStr(statistics.stdev(solver_stats[solver]["times"]))
                T.add_row([solver, successes, mean, stdev])
                T_color.add_row([solver, colors[successes != str(iterations)] + successes + reset, mean, stdev])

            print_and_log(T.get_string(), to_print=T_color.get_string(), include_brief=True)
            T.clear_rows()
            T_color.clear_rows()


if __name__ == '__main__':
    main()
