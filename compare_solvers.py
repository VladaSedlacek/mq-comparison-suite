from functools import reduce
from pathlib import Path
from subprocess import call
from prettytable import PrettyTable
import click
import datetime
import time


def secondsToStr(t):
    return "%d:%02d:%02d.%03d" % \
        reduce(lambda ll, b: divmod(ll[0], b) + ll[1:],
               [(t*1000,), 1000, 60, 60])


def print_and_log(log_paths, msg, printing=True):
    if printing:
        print(msg)
    for path in log_paths:
        with open(path, 'a') as f:
            f.write(msg + "\n")


@ click.command()
@ click.option('--o2_lb', default=6, help='lower bound for o2', type=int)
@ click.option('--o2_ub', default=20, help='upper bound for o2', type=int)
@ click.option('--runs', default=2, help='number of instances for each parameter set', type=int)
@ click.option('-v', '--verbose', default=False, is_flag=True, help='control the log verbosity')
@ click.option('-t', '--table', default=True, is_flag=True, help='display the results as a table')
def main(o2_lb, o2_ub, runs, verbose, table):
    result = ["success", "failure"]
    solvers = ['xl', 'crossbred', 'mq', 'libfes', 'wdsat', 'cms']
    q_range = [2, 16]
    o2_range = range(o2_lb, o2_ub, 2)
    log_path_1 = Path("comparison_log_brief.txt")
    log_path_2 = Path("comparison_log_verbose.txt")
    log_paths = [log_path_1, log_path_2]
    if table:
        T = PrettyTable()
        T.field_names = ["Solver", "Result", "Time"]
        T.align = "c"
    star_length = 105
    stars = '*' * star_length
    left_pad = ' ' * int((star_length - 70)/2)
    printing = not table
    print_and_log(log_paths, f"{stars}\nStarting solver comparison. See the results in {log_path_1} and {log_path_2}.")
    print_and_log(log_paths, f"Current datetime: {datetime.datetime.now().isoformat(' ', 'seconds')}")
    for seed in range(runs):
        for q in q_range:
            for o2 in o2_range:
                m = 2 * o2
                n = 3 * o2
                gen_cmd = f"sage rainbow_attacks.sage --seed {seed} --q {q} --o2 {o2} --m {m} --n {n}"
                gen_msg = f"Generating equations for seed = {seed}, q = {q}, o2 = {o2}, m = {m}, n = {n}..."
                print_and_log(log_paths, f"\n\n{stars}\n{left_pad}{gen_msg}\n{stars}")
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
                        time_seconds = secondsToStr(time_taken)
                    except Exception as e:
                        print_and_log(log_paths, str(e))
                        continue
                    # if table is True, both stdout and the brief log will be only in table-like format
                    print_and_log(log_paths[table:], f"Solver: {solver}", printing=printing)
                    print_and_log(log_paths[table:], f"Result: {result[code]}", printing=printing)
                    print_and_log(log_paths[table:], f"Time:   {time_seconds}", printing=printing)
                    print_and_log(log_paths[table:], stars, printing=printing)
                    if table:
                        T.add_row([solver, result[code], time_seconds])
                if table:
                    print_and_log(log_paths[:1], T.get_string())
                    T.clear_rows()


if __name__ == '__main__':
    main()
