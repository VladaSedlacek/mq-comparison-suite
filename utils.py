#!/usr/bin/env python3

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
