#!/usr/bin/env python3

from coastal_engineering import wave
import matplotlib.pyplot as plt
import numpy as np


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
    plt.plot(w.t, w.eta, label="Raw data")
    for m, ms in {0: "Arithmetic", 1: "Linear", 2: "Parabolic"}.items():
        plt.plot(w.t, w.eta - w.eta_correction(m), label=f"{ms} correction")
    plt.xlabel("Time (s)")
    plt.xlabel("Water level (m)")
    plt.legend()
    plt.show()


if __name__ == "__main__":
    main()

