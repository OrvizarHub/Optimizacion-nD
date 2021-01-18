import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def Main():
    data = pd.read_csv("..\debug\data.csv")

    xmax = 15

    fig, axs = plt.subplots(1, 3, figsize=(15,4))

    ax = axs[0]
    ax.axhline(0.5, color='k', linewidth=3, linestyle='--')
    ax.plot(data['cont'], data['x1'], 'bo-', linewidth=3)
    ax.set_xlabel('Iteraciones')
    ax.set_ylabel(r'$x_1$')
    ax.grid(True)
    ax.set_xlim([0, xmax])

    ax = axs[1]
    ax.axhline(0.12, color='k', linewidth=3, linestyle='--')
    ax.plot(data['cont'], data['x2'], 'bo-', linewidth=3)
    ax.set_xlabel('Iteraciones')
    ax.set_ylabel(r'$x_2$')
    ax.grid(True)
    ax.set_xlim([0, xmax])

    ax = axs[2]
    ax.axhline(-9.6547, color='k', linewidth=3, linestyle='--')
    ax.plot(data['cont'], data['cost'], 'bo-', linewidth=3)
    ax.set_xlabel('Iteraciones')
    ax.set_ylabel(r'$U$')
    ax.grid(True)
    ax.set_xlim([0, xmax])

    plt.show()


if __name__ == "__main__":
    Main()