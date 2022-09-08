from functools import reduce
from pathlib import Path
from subprocess import call
import click
import datetime
import time


def secondsToStr(t):
    return "%d:%02d:%02d.%03d" % \
        reduce(lambda ll, b: divmod(ll[0], b) + ll[1:],
               [(t*1000,), 1000, 60, 60])


def print_and_log(log_paths, msg, sep="\n"):
    print(msg)
    for path in log_paths:
        with open(path, 'a') as f:
            f.write(msg + sep)


@ click.command()
@ click.option('--o2_lb', default=6, help='lower bound for o2', type=int)
@ click.option('--o2_ub', default=20, help='upper bound for o2', type=int)
@ click.option('--runs', default=2, help='number of instances for each parameter set', type=int)
@ click.option('-v', '--verbose', default=False, is_flag=True, help='control the log verbosity')
def main(o2_lb, o2_ub, runs, verbose):
    result = ["success", "failure"]
    solvers = ['xl', 'crossbred', 'mq', 'libfes', 'wdsat', 'cms']
    q_range = [2, 16]
    o2_range = range(o2_lb, o2_ub, 2)
    log_path_1 = Path("comparison_log_brief.txt")
    log_path_2 = Path("comparison_log_verbose.txt")
    log_paths = [log_path_1, log_path_2]
    print_and_log(
        log_paths, f"Starting solver comparison, the results will be found in {log_path_1} and {log_path_2}.")
    print_and_log(
        log_paths, f"Current datetime: {datetime.datetime.now().isoformat(' ', 'seconds')}")
    for seed in range(runs):
        for q in q_range:
            for o2 in o2_range:
                m = 2 * o2
                n = 3 * o2
                generate_command = f"sage rainbow_attacks.sage --seed {seed} --q {q} --o2 {o2} --m {m} --n {n}"
                print_and_log(
                    log_paths, f"\n\n\n{'*' * 110}\nGenerating equations for seed={seed}, q={q}, o2={o2}, m={m}, n={n}...")
                call(generate_command, shell=True)
                for solver in solvers:
                    solve_command = f"sage rainbow_attacks.sage --seed {seed} --q {q} --o2 {o2} --m {m} --n {n} --solver {solver} --solve_only"
                    msg = f"\n{'*' * 110}\nExecuting: {solve_command}\n"
                    print_and_log(log_paths, msg)
                    solve_command += f" 2>> {str(log_path_2)}"
                    if verbose:
                        solve_command += f" | tee -a {str(log_path_2)}"
                    try:
                        start_time = time.time()
                        code = call(
                            solve_command + " | grep 'Attack successful!' --quiet", shell=True)
                        time_taken = time.time() - start_time
                        print_and_log(
                            log_paths, f"Result: {result[code]}")
                        print_and_log(
                            log_paths, f"Time:   {secondsToStr(time_taken)}", sep="")
                    except Exception as e:
                        print_and_log(log_paths, str(e))
                        continue


if __name__ == '__main__':
    main()
