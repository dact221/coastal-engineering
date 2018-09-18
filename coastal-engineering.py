#!/usr/bin/env python3

from coastal_engineering import wave, stats
import matplotlib.pyplot as plt
import numpy as np


METHODS = {0: "Arithmetic", 1: "Linear", 2: "Parabolic"}


REPORT = (
    "** {} **\n"
    "Mean: {:0.4f}\n"
    "Root mean square: {:0.4f}\n"
    "Skewness: {:0.4f}\n"
    "Mean of the highest half: {:0.4f}\n"
    "Mean of the highest third: {:0.4f}\n"
    "Mean of the highest tenth: {:0.4f}\n"
    "Mean of the highest hundredth: {:0.4f}\n"
    "Maximum: {:0.4f}"
)


def get_filename():
    from tkinter import Tk
    from tkinter.filedialog import askopenfilename
    root = Tk()
    root.withdraw()
    filename = askopenfilename()
    root.destroy()
    return filename


def print_report(title, arr):
    """Print statistical report.
    
    `arr` should be sorted in descending order.
    """
    print(REPORT.format(
        title,
        np.mean(arr),
        stats.rms(arr),
        stats.skw(arr),
        stats.mean2(arr),
        stats.mean3(arr),
        stats.mean10(arr),
        stats.mean100(arr),
        np.max(arr)
    ))


def main():

    # Load data
    filename = get_filename()
    print(f"Opening {filename}")
    t, eta = np.loadtxt(filename, unpack=True)
    w = wave.Wave(t, eta)

    # Plot possible signal corrections
    plt.plot(w.t, w.eta, label="Raw data")
    for m, s in METHODS.items():
        plt.plot(w.t, w.eta - w.eta_corr(m), label=f"{s} correction")
    plt.xlabel("Time (s)")
    plt.ylabel("Water level (m)")
    plt.legend()
    plt.show()

    # Perform signal correction
    m = int(input("Tidal correction method [0/1/2]: "))
    if m in [0, 1, 2]:
        w.eta -= w.eta_corr(m)

    # Characterize signal
    w.zero_crossings()

    # Plot histograms
    for i, b in enumerate([False, True]):
        plt.subplot(1, 2, i + 1)
        plt.hist(w.h, bins="doane", cumulative=b)
        plt.xlabel("Wave height (m)")
        plt.ylabel("Frequency")
    plt.show()

    # Print results
    print_report("Wave height", np.sort(w.h)[::-1])
    print_report("Period", np.sort(w.p)[::-1])


if __name__ == "__main__":
    main()

