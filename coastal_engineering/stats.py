from numpy import mean, sqrt, sum
# TODO: Make these computations more efficient

def rms(a):
    return sqrt(mean(a**2))

def skw(a):
    return mean(a**3) / rms(a)**3

def mean2(a):
    n = len(a)
    return 2 * sum(a[:n//2]) / n

def mean3(a):
    n = len(a)
    return (3 * sum(a[:n//3]) + a[n//3]) / n

def mean10(a):
    n = len(a)
    return (10 * sum(a[:n//10]) + 4 * a[n//10]) / n

def mean100(a):
    n = len(a)
    return (100 * sum(a[:n//100]) + 54 * a[n//100]) / n

