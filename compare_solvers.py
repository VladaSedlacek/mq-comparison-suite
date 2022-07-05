from subprocess import Popen
import click


@ click.command()
@ click.option('--q_lb', help='lower bound for the field size', type=int)
@ click.option('--q_ub', help='upper bound for the field size', type=int)
@ click.option('--o2_lb', default=5, help='lower bound for o2', type=int)
@ click.option('--o2_ub', default=6, help='upper bound for o2', type=int)
@ click.option('--m_lb', help='lower bound for the number of Rainbow equations', type=int)
@ click.option('--m_ub', help='upper bound for the number of Rainbow equations', type=int)
@ click.option('--n_lb', help='lower bound for the number of Rainbow variables', type=int)
@ click.option('--n_ub', help='upper bound for the number of Rainbow variables', type=int)
@ click.option('--runs', default=2, help='number of instances for each parameter set', type=int)
def main(q_lb, q_ub, o2_lb, o2_ub, m_lb, m_ub, n_lb, n_ub, runs):
    solvers = ['xl', 'crossbred', 'mq', 'libfes', 'wdsat', 'cms']
    q_range = [2, 4, 8, 16]
    o2_range = range(o2_lb, o2_ub)
    for (i, solver) in enumerate(solvers):
        for seed in range(runs):
            for q in q_range:
                for o2 in o2_range:
                    m = 2 * o2
                    n = 3 * o2
                    solve_command = f"time sage rainbow_attacks.sage --seed {seed}  --q {q} --o2 {o2} --m {m} --n {n} --solver {solver}"
                    if i > 0:
                        solve_command += " --solve_only"
                    print("\n\n" + "*" * 100)
                    print("Executing command:", solve_command)
                    Popen(solve_command, shell=True).wait()


if __name__ == '__main__':
    main()
