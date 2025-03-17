import matplotlib
matplotlib.use('Agg')
from flask import send_file
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import os.path as pt
import os




def save_plot_to_png(plot_object, filename):
    if plot_object is None:
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, 'No plot available', horizontalalignment='center', verticalalignment='center')
        plot_object = fig

    plot_filepath = os.path.join('static', 'plots', filename)
    plot_object.savefig(plot_filepath)
    plot_png = f"{filename}.png"
    return plot_png


if __name__ == '__main__':
    main()