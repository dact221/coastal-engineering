import numpy as np

def ms(a):
    return np.mean(a**2)

def rms(a):
    return np.sqrt(np.mean(a**2))

def skw(a):
    return np.mean(a**3) / rms(a)**3

def mean2(a):
    n = len(a)
    return 2 * np.sum(a[:n//2]) / n

def mean3(a):
    n = len(a)
    return (3 * np.sum(a[:n//3]) + h[n//3]) / n

def mean10(a):
    n = len(a)
    return (10 * np.sum(a[:n//10]) + 4 * h[n//10]) / n

def mean100(a):
    n = len(a)
    return (100 * np.sum(a[:n//100]) + 54 * h[n//100]) / n

def cumfreq(a, bins=3):
    hist, bin_edges = np.histogram(a, bins)
    cumhist = np.cumsum(hist)
    return cumhist, bin_edges
