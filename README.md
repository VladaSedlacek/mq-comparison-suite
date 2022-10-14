# MQ Comparison Suite

This project aims to provide a simple unified interface for state-of-the-art solvers for the multivariate quadratic ([MQ](https://eprint.iacr.org/2005/393)) problem, and to facilitate comparing their practical performance. Pull requests are welcome!

## Solvers

Currently, the following solver implementations are supported:

* [crossbred](https://eprint.iacr.org/2017/372) (the original version, not public)
* [cryptominisat](https://github.com/msoos/cryptominisat) (SAT solver for crypto problems)
* [libfes](https://github.com/cbouilla/libfes-lite) (fast exhaustive search)
* [magma](https://magma.maths.usyd.edu.au) (using the F4 algorithm for Groebner bases, closed source)
* [mq](https://gitlab.lip6.fr/almasty/mq) (simple deterministic algorithm)
* [mqsolver](https://github.com/kcning/mqsolver) (GPU version of crossbred)
* [wdsat](https://github.com/mtrimoska/WDSat) (SAT solver for systems coming from Weil descent)
* [XL](http://polycephaly.org/projects/xl) (might need a few few tweaks to run with modern compilers and Sage 9)

The solvers need to be installed externally (though none of them are required).

## Equation systems

For best results, this comparison suite should be used for equations over either $\mathbb{F}\_2$ and $\mathbb{F}\_{16}$. Note that Magma supports all finite fields, XL supports $\mathbb{F}\_2$, $\mathbb{F}\_{16}$ and $\mathbb{F}\_{31}$, while the other solvers support only $\mathbb{F}\_2$. However, using the Weil descent, we can turn any system of $m$ equations with $n$ variables over $\mathbb{F}\_{2^d}$ into a system of $dm$ equations with $dn$ variables over $\mathbb{F}\_2$. This is done by default where needed to allow comparisons.

Since some of the solvers exit after finding the first solution, this tool works for best for systems with a single solution.

While the top-level comparison is done on equation systems coming from Rainbow attacks, these can be replaced by any other systems, e.g., from the [Fukuoka challenge](https://www.mqchallenge.org/).

## Scripts
Default paths and solvers can be configured in `config_utils.py`.
The callable scripts (via click) include:
* `rainbow_attacks.sage` generates an instance of the [Rainbow](https://www.pqcrainbow.org/) cryptosystem, mounts the [differential attack](https://eprint.iacr.org/2022/214) and saves the resulting system in different formats for further usage.
  * example to just generate systems: `sage rainbow_attacks.sage --seed 0 --q 2 --o2 6 --gen_only`
  * example to also start solving: `sage rainbow_attacks.sage --seed 0 --q 2 --o2 6 --solver cms`
  * see `sage rainbow_attacks.sage --help` for more
* `invoke_solver.py` calls a solver on a provided equation file. It is called internally by `rainbow_attacks.sage`, but can also be used for custom equation systems.
  * example: `python3 invoke_solver.py --solver xl --q 2 --m 11 --n 4 -e my_equation_system`
  * see `python3 invoke_solver.py --help` for usage
* `check_solution.sage` checks the log output after a solver call and compares it with an expected solution. It is called internally by other scripts.
* `solve_fukuoka.sage` is a simple demo to wrap solving Fukuoka challenges in a uniform way.
  * example: `sage ./solve_fukuoka.sage --fukuoka_challenge_path ../fukuoka/ToyExample-type1-n10/toy_example/ToyExample-type1-n10-seed0 --solver libfes` loads the Fukuoka challenge at the given path, converts it to all other formats and then tries to solve it with libfes
* `compile_solver.py` takes care of solver which need to be precompiled in advance. It is called internally by other scripts.
* `compare_solvers.py` is the top layer wrapper which performs a comprehensive comparison of all the solvers and logs the results.
  * example: `sage compare_solvers.py --q 2 -s magma -t 3600` runs the comparison of all solvers over $\mathbb{F}\_2$ except for magma, with one hour timeout for each solver

## Equation file formats

| Solver       | Extension |
|--------------|-----------|
| `cb_gpu`     | `.cb_gpu` |
| `cb_orig`    | `.cb_orig`|
| `cms`        | `.cnf`    |
| `libfes`     | `.mq`     |
| `magma`      | `.magma`  |
| `mq`         | `.mq`     |
| `wdsat`      | `.anf`    |
| `xl`         | `.xl`     |
