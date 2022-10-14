#!/usr/bin/env python3
from math import log
from pathlib import Path


def defaults(key):
    defaults_dict = {
        "solvers": ['cb_gpu', 'cb_orig', 'cms', 'libfes', 'magma', 'mq', 'wdsat', 'xl'],
        "solvers_to_compile": ["cb_orig", "wdsat", "xl"],
        "cb_gpu_path": Path("..", "mqsolver"),
        "cb_orig_path": Path("..", "crossbred"),
        "cms_path": Path("..", "cryptominisat", "build"),
        "libfes_path": Path("..", "libfes-lite", "build"),
        "magma_path": Path("magma"),
        "mq_path": Path("..", "mq"),
        "wdsat_path": Path("..", "WDSat"),
        "xl_path": Path("..", "xl"),
        "log_path": Path("log.txt"),
        "comparison_brief": Path("comparison_log_brief.txt"),
        "comparison_verbose": Path("comparison_log_verbose.txt"),
        "timeout": 1000
    }
    return defaults_dict[key]


def solvers_to_skip():
    return ['cb_gpu', 'cb_orig', 'magma', 'xl']


def use_weil(solver_or_format):
    # determine if Weil descent should be used by default
    weil_dict = {
        'anf': True,
        'cb_gpu': True,
        'cb_orig': True,
        'cms': True,
        'cnf': True,
        'libfes': True,
        'magma': False,
        'mq': True,
        'wdsat': True,
        'xl': False,
    }
    return weil_dict[solver_or_format]


def get_eq_format(solver):
    formats = {
        'cb_gpu': 'cb_gpu',
        'cb_orig': 'cb_orig',
        'cms': 'cnf',
        'libfes': 'mq',
        'magma': 'magma',
        'mq': 'mq',
        'wdsat': 'anf',
        'xl': 'xl',
    }
    return formats[solver]


def get_log_format(solver):
    formats = {
        'cb_gpu': 'cb',
        'cb_orig': 'cb',
        'cms': 'cms',
        'libfes': 'mq',
        'magma': 'magma',
        'mq': 'mq',
        'wdsat': 'wdsat',
        'xl': 'xl',
    }
    return formats[solver]


def get_rainbow_dims(o2, m=None, n=None):
    if m is None:
        m = 2 * o2
    if n is None:
        n = 3 * o2
    M = m - 1
    N = n - m - 2
    return m, n, M, N


def get_weil_dims(q, M, N):
    assert q % 2 == 0
    d = int(log(q, 2))
    return M * d, N * d


def get_dim_str(q, M, N, weil):
    assert q % 2 == 0
    d = int(log(q, 2))
    return f"_M_{M * d}_N_{N * d}_weil" if (d > 1 and weil) else f"_M_{M}_N_{N}"


def declare_paths(seed, q, o2, m, n):
    system_folder_path = Path("systems")
    Path(system_folder_path).mkdir(parents=True, exist_ok=True)
    base_system_name = f"rainbow_diff_s_{seed}_q_{q}_o2_{o2}_m_{m}_n_{n}"
    return system_folder_path, base_system_name


def get_sol_path(seed, q, o2, m, n, M, N, weil):
    dim_str = get_dim_str(q, M, N, weil)
    system_folder_path, base_system_name = declare_paths(seed, q, o2, m, n)
    return Path(system_folder_path, f"{base_system_name}{dim_str}.sol")


def get_eq_path(seed, q, o2, m, n, M, N, solver):
    dim_str = get_dim_str(q, M, N, use_weil(solver))
    system_folder_path, base_system_name = declare_paths(seed, q, o2, m, n)
    return Path(system_folder_path, f"{base_system_name}{dim_str}.{get_eq_format(solver)}")
