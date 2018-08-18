from numpy import loadtxt
import waves

_, t, eta = loadtxt(
    "data2.csv",
    delimiter=",",
    skiprows=1,
    unpack=True
)

w = waves.Wave(t, eta)
w.zero_crossings(0)

print(f"Periods: {w.p}")
print(f"Wave heights: {w.h}")
