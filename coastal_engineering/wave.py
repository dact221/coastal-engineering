import numpy as np

class Wave:

    def __init__(self, t, eta, verbose=False):
        self.t = t
        self.eta = eta
        self.verbose = verbose
        self.n = len(eta)

    def eta_correction(self, method):
        if method == 0:
            av = np.mean(self.eta)
        elif method == 1:
            av = self._lin()
        elif method == 2:
            av = self._parab()
        self.eta -= av

    def _lin(self):
        ns = np.arange(self.n) + 1
        n0, n1, n2 = self._n(2, ns)
        y0, y1 = self._y(1, ns)
        a0 = (n2 * y0 - n1 * y1) / (n0 * n2 - n1**2)
        a1 = (n0 * y1 - n1 * y0) / (n0 * n2 - n1**2)
        return a0 + a1 * ns

    def _parab(self):
        ns = np.arange(self.n) + 1
        n0, n1, n2, n3, n4 = self._n(4, ns)
        y0, y1, y2 = self._y(2, ns)
        d = n0*n2*n4 + 2*n1*n2*n3 - n2**3 - n0*n3**2 - n4*n1**2
        b0 = (y0*(n2*n4 - n3**2) + y1*(n2*n3 - n1*n4) + y2*(n1*n3 - n2**2)) / d
        b1 = (y0*(n2*n3 - n1*n4) + y1*(n0*n4 - n2**2) + y2*(n1*n2 - n0*n3)) / d
        b2 = (y0*(n1*n3 - n2**2) + y1*(n1*n2 - n0*n3) + y2*(n0*n2 - n1**2)) / d
        return b0 + b1*ns + b2*ns**2

    def _n(self, r, ns):
        exp = np.arange(r + 1)[:, np.newaxis]
        return np.sum(ns**exp, axis=1)

    def _y(self, r, ns):
        exp = np.arange(r + 1)[:, np.newaxis]
        return np.sum(self.eta * ns**exp, axis=1)

    def _roots(self, iz):
        m = (self.t[iz+1] - self.t[iz]) / (self.eta[iz+1] - self.eta[iz])
        return -self.eta[iz] * m + self.t[iz]

    def _max_parab(self, i):
        a = (self.eta[i-1] - 2 * self.eta[i] + self.eta[i+1]) / 2
        b = (self.eta[i+1] - self.eta[i-1]) / 2
        c = self.eta[i]
        return c - b**2 / 4 / a

    def _zero(self, method):
        a = self.eta[:-1]
        b = self.eta[1:]
        if method == 0:
            return np.where((a < 0) & (b > 0))[0]
        elif method == 1:
            return np.where((a > 0) & (b < 0))[0]

    def zero_crossings(self, method):
        iz = self._zero(method)
        r = self._roots(iz)

        self.p = r[1:] - r[:-1] # periods
        # TODO: Improve extrema finding
        iext = np.empty((2, len(self.p)), dtype=int)
        for i, (a, b) in enumerate(zip(iz[:-1] + 1, iz[1:])):
            iext[0][i] = np.argmin(self.eta[a:b]) + a
            iext[1][i] =  np.argmax(self.eta[a:b]) + a
        eta_min, eta_max = self._max_parab(iext)

        self.h = eta_max - eta_min # wave heights

        if self.verbose == True:
            print("Zero crossings:", iz)
            print("Roots:", r)
            print("i minimum:", i[0])
            print("i maximum:", i[1])
            print("eta min.:", eta_min)
            print("eta max.:", eta_max)

