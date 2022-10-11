#!/usr/bin/env python3

import matplotlib.pyplot as plt
from pathlib import Path
import json


json_path = Path("comparison.json")
with open(json_path) as j:
    results = json.load(j)


solvers = ['cb_gpu', 'cb_orig', 'cms', 'libfes', 'magma', 'mq', 'wdsat', 'xl']
q = 2
iterations = 2
max_suc = f"{iterations} of {iterations}"
# pick only those o2 for which the computations have finished
o2s = [o2 for o2 in results[f"q={q}"].keys() if any(results[f"q={q}"][o2].values())]
o2x = [int(s.split("=")[1]) for s in o2s]


def stats(o2, solver, q=2):
    return results[f"q={q}"][o2][solver]


fig = plt.figure()
for solver in solvers:
    try:
        finished = [st for st in [stats(o2, solver) for o2 in o2s] if st is not None]
    except KeyError:
        continue
    times = [par["mean_time"] if par["successes"] == max_suc else -100 for par in finished]
    plt.plot(o2x, times, label=solver)

plt.legend(loc='upper left')
plt.show()
max_suc
