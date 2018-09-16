#!/usr/bin/env python3

from coastal_engineering import wave, stats
import matplotlib.pyplot as plt
import numpy as np


REPORT = (
    "** {0} **\n"
    "Mean: {1:0.4f}\n"
    "Root mean square: {2:0.4f}\n"
    "Skewness: {3:0.4f}\n"
    "Mean of the highest half: {4:0.4f}\n"
    "Mean of the highest third: {5:0.4f}\n"
    "Mean of the highest tenth: {6:0.4f}\n"
    "Mean of the highest hundredth: {7:0.4f}\n"
    "Maximum: {8:0.4f}"
)


def get_filename():
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename
    root = Tk()
    root.withdraw()
    filename = askopenfilename()
    root.destroy()
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
    plt.ylabel("Water level (m)")
    plt.legend()
    plt.show()

    m = input("Tidal correction method [0/1/2]: ")
    w.eta -= w.eta_correction(int(m))
    w.zero_crossings()
    plt.hist(w.h, bins="doane")
    plt.show()

    it = {
        "Wave height": np.flipud(np.sort(w.h)),
        "Period": np.flipud(np.sort(w.p))
    }
    for s, arr in it.items():
        print(REPORT.format(
            s,
            np.mean(arr),
            stats.rms(arr),
            stats.skw(arr),
            stats.mean2(arr),
            stats.mean3(arr),
            stats.mean10(arr),
            stats.mean100(arr),
            np.max(arr)
        ))

if __name__ == "__main__":
    main()

