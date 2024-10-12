# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt

def plot_fit(x_data, y_data, x_fit, y_fit, ):
    
    fig, ax = plt.subplots()
    
    
    ax.plot(x_fit, y_fit, color='blue')
    ax.scatter(x_data, y_data, color='red')
    
    return fig, ax


