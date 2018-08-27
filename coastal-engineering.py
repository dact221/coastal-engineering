#!/usr/bin/env python3

import numpy as np
from coastal_engineering import wave


def get_filename():
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename
    Tk().withdraw()
    filename = askopenfilename()
    return filename


def main():
    filename = get_filename()
    print(f"Opening {filename}")
    t, eta = np.loadtxt(filename, unpack=True)

    w = wave.Wave(t, eta)
    w.zero_crossings(0) # rising zero-crossings

    print(
        f"Periods: {w.p}",
        f"Wave heights: {w.h}",
        f"Mean sea level {np.mean(w.eta)}",
        f"Sea level variance {np.var(w.eta)}",
        sep="\n"
    )


if __name__ == "__main__":
    main()
