from functools import reduce
from operator import itemgetter
from pathlib import Path
from subprocess import call
from prettytable import PrettyTable
import click
import datetime
import statistics
import time


def secondsToStr(t):
    return "%d:%02d:%02d.%03d" % \
        reduce(lambda ll, b: divmod(ll[0], b) + ll[1:],
               [(t*1000,), 1000, 60, 60])


def print_and_log(log_paths, to_log, to_print=None):
    if to_print == None:
        to_print = to_log
    if not to_print == "":
        print(to_print)
    for path in log_paths:
        with open(path, 'a') as f:
            f.write(to_log + "\n")


@ click.command()
@ click.option('--o2_lb', default=6, help='lower bound for o2', type=int)
@ click.option('--o2_ub', default=20, help='upper bound for o2', type=int)
@ click.option('--runs', default=2, help='number of instances for each parameter set', type=int)
@ click.option('-v', '--verbose', default=False, is_flag=True, help='control the log verbosity')
@ click.option('-t', '--table', default=True, is_flag=True, help='display the results as a table')
def main(o2_lb, o2_ub, runs, verbose, table):

    # Colors (from https://stackoverflow.com/questions/43583847/python-pretty-table-with-color-output)
    red = "\033[0;31;40m"
    green = "\033[0;32;40m"
    reset = "\033[0m"
    colors = [green, red]

    results = ["success", "failure"]
    solvers = ['xl', 'crossbred', 'mq', 'libfes', 'wdsat', 'cms']
    q_range = [2, 16]
    o2_range = range(o2_lb, o2_ub, 2)
    log_path_1 = Path("comparison_log_brief.txt")
    log_path_2 = Path("comparison_log_verbose.txt")
    log_paths = [log_path_1, log_path_2]
    if table:
        T = PrettyTable()
        res_col_name = f"Successes (out of {runs})"
        T.field_names = ["Solver", f"{res_col_name}", "Average time", "Standard deviation of time"]
        T.align = "c"

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

    star_length = 105
    stars = '*' * star_length
    left_pad = ' ' * int((star_length - 70)/2)
    print_and_log(log_paths, f"{stars}\nStarting solver comparison. See the results in {log_path_1} and {log_path_2}.")
    print_and_log(log_paths, f"Current datetime: {datetime.datetime.now().isoformat(' ', 'seconds')}")
    for q in q_range:
        for o2 in o2_range:
            m = 2 * o2
            n = 3 * o2
            gen_msg = f"Generating equations for q = {q}, o2 = {o2}, m = {m}, n = {n}; {runs} iterations..."
            print_and_log(log_paths, f"\n\n{stars}\n{left_pad}{gen_msg}\n{stars}")
            solver_stats = {solver: {"successes": 0, "times": []} for solver in solvers}
            for seed in range(runs):
                gen_cmd = f"sage rainbow_attacks.sage --seed {seed} --q {q} --o2 {o2} --m {m} --n {n}"
                call(gen_cmd, shell=True)
                for solver in solvers:
                    solve_cmd = f"sage rainbow_attacks.sage --seed {seed} --q {q} --o2 {o2} --m {m} --n {n} --solver {solver} --solve_only"
                    if verbose:
                        solver_verbose = f"\n{stars}\nExecuting: {solve_cmd}\n"
                        print_and_log([log_path_2], solver_verbose)
                        solve_cmd += f" 2 >> {str(log_path_2)} | tee -a {str(log_path_2)}"
                    else:
                        solve_cmd += f" 2>> {str(log_path_2)}"
                    try:
                        start_time = time.time()
                        code = call(solve_cmd + " | grep 'Attack successful!' --quiet", shell=True)
                        time_taken = time.time() - start_time
                    except Exception as e:
                        print_and_log(log_paths, str(e))
                        continue
                    solver_stats[solver]["successes"] += (1-code)
                    solver_stats[solver]["times"].append(time_taken)

                    # if table is True, both stdout and the brief log will be only in table-like format
                    print_and_log(log_paths[table:], f"Solver: {solver}", to_print="")
                    print_and_log(log_paths[table:], f"Result: {results[code]}", to_print="")
                    print_and_log(log_paths[table:], f"Time:   {secondsToStr(time_taken)}", to_print="")
                    print_and_log(log_paths[table:], stars, to_print="")

            for solver in solvers:
                successes = str(solver_stats[solver]["successes"])
                mean = secondsToStr(statistics.mean(solver_stats[solver]["times"]))
                stdev = secondsToStr(statistics.stdev(solver_stats[solver]["times"]))

                if table:
                    T.add_row([solver, successes, mean, stdev])
                    T_color.add_row([solver, colors[successes != str(runs)] + successes + reset, mean, stdev])
            if table:
                print_and_log(log_paths[:1], T.get_string(), to_print=T_color.get_string())
                T.clear_rows()
                T_color.clear_rows()


if __name__ == '__main__':
    main()
