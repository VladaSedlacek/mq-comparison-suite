# MQ Comparison Suite

This project aims to provide a simple unified interface for state-of-the-art solvers for the multivariate quadratic ([MQ](https://eprint.iacr.org/2005/393)) problem, and to facilitate comparing their practical performance.

## Solvers

Currently, the following solver implementations are supported:

* [crossbred](https://eprint.iacr.org/2017/372) (the original version, not public)
* [cryptominisat](https://github.com/msoos/cryptominisat) (SAT solver for crypto problems)
* [libfes](https://github.com/cbouilla/libfes-lite) (fast exhaustive search)
* [magma](https://magma.maths.usyd.edu.au) (using the F4 algorithm for Groebner bases, closed source)
* [mq](https://gitlab.lip6.fr/almasty/mq) (simple deterministic algorithm)
* [mqsolver](https://github.com/kcning/mqsolver) (GPU version of crossbred)
* [wdsat](https://github.com/mtrimoska/WDSat) (SAT solver for systems coming from Weil descent)
* [XL](http://polycephaly.org/projects/xl)

The solvers need to be installed externally (though none of them are required).

Note that magma supports all finite fields, XL supports $\mathbb{F}\_2$, $\mathbb{F}\_{16}$ and $\mathbb{F}\_{31}$, while the other solvers support only $\mathbb{F}\_2$. However, using the Weil descent, we can turn any system of $m$ equations with $n$ variables over $\mathbb{F}\_{2^r}$ into a system of $rm$ equations with $rn$ variables over $\mathbb{F}\_2$. This is done by default where needed to allow comparisons.

## Scripts

* `rainbow_attacks.sage` generates an instance of the [Rainbow](https://www.pqcrainbow.org/) cryptosystem, mounts the [differential attack](https://eprint.iacr.org/2022/214) and saves the resulting system in different formats for further usage.
  * example to just generate systems: `sage rainbow_attacks.sage --seed 0 --q 2 --o2 6 --gen_only`
  * example to also start solving: `sage rainbow_attacks.sage --seed 0 --q 2 --o2 6 --solver cms`
  * see `sage rainbow_attacks.sage --help` for more
* `invoke_solver.py` calls a solver on a provided equation file. It is called internally by `rainbow_attacks.sage`, but can also be used for custom equation systems.
  * example: `python3 invoke_solver.py --solver xl --q 2 --m 11 --n 4 -e my_equation_system`
  * see `python3 invoke_solver.py --help` for usage
* `compile_solver.py` takes care of solver which need to be precompiled in advance. It is called internally by other scripts.
* `compare_solvers.py` is the top layer wrapper which performs a comprehensive comparison of all the solvers and logs the results.
  * example: `sage compare_solvers.py -s magma` runs the comparison of all solvers except for magma

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
