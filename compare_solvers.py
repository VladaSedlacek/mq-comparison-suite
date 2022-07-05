from pathlib import Path
from subprocess import call
import click


def print_and_log(f, msg):
    print(msg)
    f.write(msg + "\n\n")


@ click.command()
@ click.option('--o2_lb', default=5, help='lower bound for o2', type=int)
@ click.option('--o2_ub', default=6, help='upper bound for o2', type=int)
@ click.option('--runs', default=2, help='number of instances for each parameter set', type=int)
@ click.option('-v', '--verbose', default=False, is_flag=True, help='control the output verbosity')
def main(o2_lb, o2_ub, runs, verbose):
    result = ["success", "failure"]
    solvers = ['xl', 'crossbred', 'mq', 'libfes', 'wdsat', 'cms']
    q_range = [2, 4, 8, 16]
    o2_range = range(o2_lb, o2_ub)
    log_path = Path("comparison_log.txt")
    with open(log_path, 'w') as file:
        for seed in range(runs):
            for q in q_range:
                for o2 in o2_range:
                    m = 2 * o2
                    n = 3 * o2
                    for (i, solver) in enumerate(solvers):
                        solve_command = f"sage rainbow_attacks.sage --seed {seed} --q {q} --o2 {o2} --m {m} --n {n} --solver {solver}"
                        if i > 0:
                            solve_command += " --solve_only"
                        if verbose:
                            solve_command += " | tee -a " + str(log_path)
                        msg = f"\n\n{'*' * 100}\nExecuting command: {solve_command}\n"
                        print_and_log(file, msg)
                        try:
                            code = call(solve_command +
                                        " | grep 'Attack successful!'", shell=True)
                            print_and_log(file,  f"RESULT: {result[code]}")
                        except Exception as e:
                            print_and_log(file, str(e))
                            continue


if __name__ == '__main__':
    main()
