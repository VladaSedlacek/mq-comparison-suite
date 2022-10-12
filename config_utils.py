#!/usr/bin/env python3
from math import log
from pathlib import Path


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


def declare_paths(seed, q, o2, m, n):
    system_folder_path = Path("systems")
    Path(system_folder_path).mkdir(parents=True, exist_ok=True)
    base_system_name = f"rainbow_diff_s_{seed}_q_{q}_o2_{o2}_m_{m}_n_{n}"
    return system_folder_path, base_system_name


def get_sol_path(seed, q, o2, m, n, M, N, weil=True):
    assert q % 2 == 0
    d = int(log(q, 2))
    dim_str = f"_M_{M * d}_N_{N * d}_weil" if (d > 1 and weil) else f"_M_{M}_N_{N}"
    system_folder_path, base_system_name = declare_paths(seed, q, o2, m, n)
    return Path(system_folder_path, f"{base_system_name}{dim_str}.sol")


def get_eq_path(seed, q, o2, m, n, M, N, solver):
    assert q % 2 == 0
    d = int(log(q, 2))
    dim_str = f"_M_{M * d}_N_{N * d}_weil" if (d > 1 and use_weil(solver)) else f"_M_{M}_N_{N}"
    system_folder_path, base_system_name = declare_paths(seed, q, o2, m, n)
    return Path(system_folder_path, f"{base_system_name}{dim_str}.{get_eq_format(solver)}")
