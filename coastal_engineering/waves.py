import numpy as np

class Wave:

    def __init__(self, t, eta, verbose=False):
        self.t = t
        self.eta = eta
        self.verbose = verbose

    def mean_lvl_correction(self, m):
        if m == 0:
            av = np.mean(self.water_lvl)
        elif m == 1:
            av = self.lin()
        elif m == 2:
            av = self.parab()
        self.eta -= av

    def lin(self):
        ns = np.arange(len(self.water_lvl)) + 1
        n0, n1, n2 = self.n(2, ns)
        y0, y1 = self.y(1, ns)
        a0 = (n2 * y0 - n1 * y1) / (n0 * n2 - n1**2)
        a1 = (n0 * y1 - n1 * y0) / (n0 * n2 - n1**2)
        return a0 + a1 * ns

    def parab(self):
        ns = np.arange(len(self.water_lvl)) + 1
        n0, n1, n2, n3, n4 = self.n(4, ns)
        y0, y1, y2 = self.y(2, ns)
        d = n0*n2*n4 + 2*n1*n2*n3 - n2**3 - n0*n3**2 - n4*n1**2
        b0 = (y0*(n2*n4 - n3**2) + y1*(n2*n3 - n1*n4) + y2*(n1*n3 - n2**2)) / d
        b1 = (y0*(n2*n3 - n1*n4) + y1*(n0*n4 - n2**2) + y2*(n1*n2 - n0*n3)) / d
        b2 = (y0*(n1*n3 - n2**2) + y1*(n1*n2 - n0*n3) + y2*(n0*n2 - n1**2)) / d
        return b0 + b1*ns + b2*ns**2

    def n(self, r, ns):
        exp = np.arange(r + 1).reshape((-1, 1))
        return np.sum(ns**exp, axis=1)

    def y(self, r, ns):
        exp = np.arange(r + 1).reshape((-1, 1))
        return np.sum(self.eta * ns**exp, axis=1)

    def roots(self, iz):
        m = (self.t[iz+1] - self.t[iz]) / (self.eta[iz+1] - self.eta[iz])
        return -self.eta[iz] * m + self.t[iz]

    def max_parab(self, i):
        a = (self.eta[i-1] - 2 * self.eta[i] + self.eta[i+1]) / 2
        b = (self.eta[i+1] - self.eta[i-1]) / 2
        c = self.eta[i]
        return c - b**2 / 4 / a

    def zero_crossings(self, m):
        if m == 0: # rising
            iz = np.where((self.eta[:-1] < 0) & (self.eta[1:] > 0))[0]
        elif m == 1: # falling
            iz = np.where((self.eta[:-1] > 0) & (self.eta[1:] < 0))[0]
        r = self.roots(iz)
        self.p = r[1:] - r[:-1] # periods
        i = np.empty((2, len(self.p)), dtype=int)
        for j, (a, b) in enumerate(zip(iz[:-1] + 1, iz[1:])):
            i[0, j] = np.argmin(self.eta[a:b]) + a
            i[1, j] =  np.argmax(self.eta[a:b]) + a
        eta_min, eta_max = self.max_parab(i)
        self.h = eta_max - eta_min # wave heights
        if self.verbose == True:
            print("Zero crossings:", iz)
            print("Roots:", r)
            print("i minimum:", i[0])
            print("i maximum:", i[1])
            print("eta min.:", eta_min)
            print("eta max.:", eta_max)

    def rising_zero_crossings_legacy(self):
        print("Using legacy version!")
        has_zero = (self.eta[:-1] < 0) & (self.eta[1:] > 0)
        rts = list()
        izero = list()
        for i, x in enumerate(has_zero):
            if x:
                m = (self.t[i+1] - self.t[i]) / (self.eta[i+1] - self.eta[i])
                r = -self.eta[i] * m + self.t[i]
                rts.append(r)
                izero.append(i)
        rts = np.array(rts)
        imax = list()
        imin = list()
        d = izero[0] + 1
        for a, b in zip(izero[:-1], izero[1:]):
            a += 1
            imax.append(np.argmax(self.eta[a:b]) + d)
            imin.append(np.argmin(self.eta[a:b]) + d)
            d = b + 1
        maxim = list()
        minim = list()
        for i in imax:
            a = (self.eta[i-1] - 2 * self.eta[i] + self.eta[i+1]) / 2
            b = (self.eta[i+1] - self.eta[i-1]) / 2
            c = self.eta[i]
            eta_max = c - b**2 / 4 / a
            t_max = self.t[i] - self.t[1] * b / 2 / a
            maxim.append((t_max, eta_max))
        for i in imin:
            a = (self.eta[i-1] - 2 * self.eta[i] + self.eta[i+1]) / 2
            b = (self.eta[i+1] - self.eta[i-1]) / 2
            c = self.eta[i]
            eta_min = c - b**2 / 4 / a
            t_min = self.t[i] - self.t[1] * b / 2 / a
            minim.append((t_min, eta_min))
        maxim = np.array(maxim)
        minim = np.array(minim)
        self.p = rts[1:] - rts[:-1]
        self.h = maxim[:, 1] - minim[:, 1]
