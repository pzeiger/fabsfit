# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np

def plot_fit(xdata, ydata, fitfunc ):
    
    fig, ax = plt.subplots()
    
    
    s = np.linspace(0, 6, 601)
    ax.plot(s, fitfunc(s), color='blue')
    ax.scatter(xdata, ydata, color='red', marker='x')
    
    return fig, ax

