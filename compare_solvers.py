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


@ click.command()
@ click.option('--o2_lb', default=6, help='lower bound for o2', type=int)
@ click.option('--o2_ub', default=20, help='upper bound for o2', type=int)
@ click.option('--runs', default=2, help='number of instances for each parameter set', type=int)
@ click.option('--log_path_brief', default=Path("comparison_log_brief.txt"), help='the path to the brief log')
@ click.option('--log_path_verbose', default=Path("comparison_log_verbose.txt"), help='the path to the verbose log')
def main(o2_lb, o2_ub, runs, log_path_brief, log_path_verbose):

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

    # Colors (from https://stackoverflow.com/questions/43583847/python-pretty-table-with-color-output)
    red = "\033[0;31;40m"
    green = "\033[0;32;40m"
    reset = "\033[0m"
    colors = [green, red]

    results = ["success", "failure"]
    solvers = ['xl', 'crossbred', 'mq', 'libfes', 'wdsat', 'cms']
    # solvers = ['mq', 'libfes', 'wdsat', 'cms']
    q_range = [2, 16]
    o2_range = range(o2_lb, o2_ub, 2)
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
    print_and_log(
        f"{stars}\nStarting solver comparison, check {log_path_verbose} for the detailed results.", include_brief=True)
    print_and_log(f"Current datetime: {datetime.datetime.now().isoformat(' ', 'seconds')}", include_brief=True)
    for q in q_range:
        for o2 in o2_range:
            m = 2 * o2
            n = 3 * o2
            gen_msg = f"Generating equations for q = {q}, o2 = {o2}, m = {m}, n = {n}; {runs} iterations..."
            print_and_log(f"\n\n{stars}\n{left_pad}{gen_msg}\n{stars}", include_brief=True)
            solver_stats = {solver: {"successes": 0, "times": []} for solver in solvers}
            for seed in range(runs):
                gen_cmd = f"sage rainbow_attacks.sage --seed {seed} --q {q} --o2 {o2} --m {m} --n {n}"
                call(gen_cmd, shell=True)
                for solver in solvers:
                    solve_cmd = f"sage rainbow_attacks.sage --seed {seed} --q {q} --o2 {o2} --m {m} --n {n} --solver {solver} --solve_only"
                    print_and_log(f"\n{stars}\nExecuting: {solve_cmd}\n", to_print="")
                    solve_cmd += f" 2>> {log_path_verbose} | tee -a {str(log_path_verbose)}"
                    try:
                        start_time = time.time()
                        code = call(solve_cmd + " | grep 'Attack successful!' --quiet", shell=True)
                        time_taken = time.time() - start_time
                    except Exception as e:
                        print_and_log(str(e))
                        continue
                    solver_stats[solver]["successes"] += (1-code)
                    solver_stats[solver]["times"].append(time_taken)

                    print_and_log(f"Solver: {solver}", to_print="")
                    print_and_log(f"Result: {results[code]}", to_print="")
                    print_and_log(f"Time:   {secondsToStr(time_taken)}", to_print="")
                    print_and_log(stars, to_print="")

            for solver in solvers:
                successes = str(solver_stats[solver]["successes"])
                mean = secondsToStr(statistics.mean(solver_stats[solver]["times"]))
                stdev = secondsToStr(statistics.stdev(solver_stats[solver]["times"]))

                T.add_row([solver, successes, mean, stdev])
                T_color.add_row([solver, colors[successes != str(runs)] + successes + reset, mean, stdev])

            print_and_log(T.get_string(), to_print=T_color.get_string(), include_brief=True)
            T.clear_rows()
            T_color.clear_rows()


if __name__ == '__main__':
    main()
