import numpy as np

class Wave:
    """This class represents a wind wave.

    In order to create a `Wave` instance you need to provide time and water
    level data as NumPy arrays. It implements methods to perform tidal
    correction and signal characterization, and to estimate parameters for the
    statistical distribution of wind waves.
    
    `p` and `h` attributes are computed by `zero_crossings` (initially `None`).

    Attributes:
        t: NumPy array containing time series.
        eta: NumPy array containing water level series.
        p: NumPy array containing spatial periods.
        h: NumPy array containing wave heights.
        n: number of data entries.
    """


    def __init__(self, t, eta):
        """Inits `Wave` with time and water level arrays"""
        self.t = t
        self.eta = eta
        self.p = None
        self.h = None
        self.n = len(eta)


    def _n(self, r):
        n = list()
        for j in range(r + 1):
            s = 0
            for i in range(1, self.n + 1):
                s += i ** j
            n.append(s)
        return n


    def _y(self, r):
        y = list()
        for j in range(r + 1):
            s = 0
            for i in range(1, self.n + 1):
                s += self.eta[i - 1] * i ** j
            y.append(s)
        return y


    def _lin(self):
        ns = np.arange(self.n) + 1
        n0, n1, n2 = self._n(2)
        y0, y1 = self._y(1)
        a0 = (n2 * y0 - n1 * y1) / (n0 * n2 - n1**2)
        a1 = (n0 * y1 - n1 * y0) / (n0 * n2 - n1**2)
        return a0 + a1 * ns


    def _parab(self):
        ns = np.arange(self.n) + 1
        n0, n1, n2, n3, n4 = self._n(4)
        y0, y1, y2 = self._y(2)
        d = n0*n2*n4 + 2*n1*n2*n3 - n2**3 - n0*n3**2 - n4*n1**2
        b0 = (y0*(n2*n4 - n3**2) + y1*(n2*n3 - n1*n4) + y2*(n1*n3 - n2**2)) / d
        b1 = (y0*(n2*n3 - n1*n4) + y1*(n0*n4 - n2**2) + y2*(n1*n2 - n0*n3)) / d
        b2 = (y0*(n1*n3 - n2**2) + y1*(n1*n2 - n0*n3) + y2*(n0*n2 - n1**2)) / d
        return b0 + b1*ns + b2*ns**2


    def eta_corr(self, method):
        """Computes tidal correction for water level data.

        Args:
            method: integer corresponding to the method to compute tidal
            correction. 0 (arithmetic), 1 (linear) or 2 (parabolic).

        Returns:
            A NumPy array containing the tidal corrections to be subtracted
            from water level entries.
        """
        if method == 0:
            av = np.mean(self.eta)
        elif method == 1:
            av = self._lin()
        elif method == 2:
            av = self._parab()
        return av


    def _zero(self, method):
        """Finds rising or falling zero-crossings in water level array.

        Args:
            method: integer. 0 for rising zero-crossings or 1 for falling
            zero-crossings.

        Returns:
            A NumPy array containing the indexes where zero-crossings occur
            in water level array.
        """
        a = self.eta[:-1]
        b = self.eta[1:]
        # TODO: What to do with zero-crossings where `eta[i] * eta[i+1] == 0`?
        if method == 0:
            return np.where((a < 0) & (b > 0))[0]
        elif method == 1:
            return np.where((a > 0) & (b < 0))[0]


    def _roots(self, i):
        """Compute roots of the wave using linear interpolation.

        Args:
            i: array of indexes where zero-crossings occur.

        Returns:
            An array containing the times at which zero-crossings occur.
        """
        m = (self.t[i+1] - self.t[i]) / (self.eta[i+1] - self.eta[i])
        return -self.eta[i] * m + self.t[i]


    def _max_parab(self, i):
        """Computes the vertex of the adjusted parabolas at the maximum and
        minimum water levels of each single wave.

        Args:
            i: a 2-row NumPy array containing the indexes where minimum and
            maximum water levels are for each single wave.

        Returns:
            A 2-row array containing adjusted minimum/maximum water levels.
        """
        a = (self.eta[i-1] - 2 * self.eta[i] + self.eta[i+1]) / 2
        b = (self.eta[i+1] - self.eta[i-1]) / 2
        c = self.eta[i]
        return c - b**2 / 4 / a


    def zero_crossings(self, method=0):
        """Computes/sets spatial periods and wave heights using rising or
        falling zero-crossings method.

        Creates or overwrite `p` (spatial periods) and `h` (wave heights)
        attributes for `Wave`.

        Args:
            method: integer. 0 for rising zero-crossings or 1 for falling
            zero-crossings.
        """
        iz = self._zero(method)
        r = self._roots(iz)

        self.p = r[1:] - r[:-1] # periods

        iext = np.empty((2, len(self.p)), dtype=int)
        for i, (a, b) in enumerate(np.column_stack((iz[:-1], iz[1:])) + 1):
            iext[0, i] = np.argmin(self.eta[a:b]) + a
            iext[1, i] =  np.argmax(self.eta[a:b]) + a
        eta_min, eta_max = self._max_parab(iext)

        self.h = eta_max - eta_min # wave heights

