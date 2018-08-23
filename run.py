import numpy as np
import waves

_, t, eta = np.loadtxt(
    "data2.csv",
    delimiter=",",
    skiprows=1,
    unpack=True
)

w = waves.Wave(t, eta)
w.zero_crossings(0)

print(f"Periods: {w.p}")
print(f"Wave heights: {w.h}")
print("Mean sea level", np.mean(w.eta))
print("Variance", np.var(w.eta))
print("Maximum", np.max(w.eta))
