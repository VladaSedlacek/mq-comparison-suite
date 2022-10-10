#!/usr/bin/env python3
from math import log
from pathlib import Path


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


def declare_paths(seed, q, o2, m, n):
    system_folder_path = Path("systems")
    Path(system_folder_path).mkdir(parents=True, exist_ok=True)
    base_system_name = f"rainbow_diff_s_{seed}_q_{q}_o2_{o2}_m_{m}_n_{n}"
    return system_folder_path, base_system_name


def get_eq_pathname(seed, q, o2, m, n, M, N, solver):
    assert q % 2 == 0
    d = int(log(q, 2))
    dim_str = f"_M_{M * d}_N_{N * d}_weil" if d > 1 else f"_M_{M}_N_{N}"
    system_folder_path, base_system_name = declare_paths(seed, q, o2, m, n)
    return Path(system_folder_path, f"{base_system_name}{dim_str}.{get_eq_format(solver)}")
