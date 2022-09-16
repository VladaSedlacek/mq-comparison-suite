#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import json


json_path = Path("comparison.json")
with open(json_path) as j:
    results = json.load(j)


solvers = ['cb_gpu', 'cb_orig', 'cms', 'libfes', 'magma', 'mq', 'wdsat', 'xl']
q = 2
iterations = 2
max_suc = f"{iterations} of {iterations}"
o2s = results[f"q={q}"].keys()
o2x = [int(s.split("=")[1]) for s in o2s]


def stats(o2, solver, q=2):
    return results[f"q={q}"][o2][solver]


fig = plt.figure()
for solver in solvers:
    times = [stats(o2, solver)["mean_time"] if stats(o2, solver)["successes"] == max_suc else 0 for o2 in o2s]
    plt.plot(o2x, times, label=solver)

plt.legend(loc='upper left')
plt.show()
max_suc
