#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
# Created by AndresD at 24/05/19

Features:
    + Plotting library for GRASS-PYTHON tools

@author:    Andres Felipe Duque Perez
Email:      aduquep@ingenevo.com.co
"""


# %% plot timeseries
def plt_timeseries(timeseries, title, var):
    """Plot pandas timeseries:\n

       title        : plot title\n
       timeseries   : pandas timeseries

    """
    # import required modules
    import matplotlib.pyplot as plt

    # new figure
    fig = plt.figure(figsize=[10, 8])
    ax = fig.add_subplot(111)

    # set axis
    ax.set_xlabel('Fecha', fontsize=12, fontweight='bold')
    if var == 1:
        ax.set_ylabel('Precipitación (mm)',
                      fontsize=12, fontweight='bold')
    else:
        ax.set_ylabel('Caudal (' + '$\mathregular{m^3}$' + '/s)',
                      fontsize=12, fontweight='bold')

    # plt timeseries
    timeseries.plot()
    plt.title(title, fontweight='bold')

    return fig


# %% plot annual cycle boxplot
def plt_annual_cycle(annual_cycle_df, annual_cycle_error_df, title, var):
    """Plot annual cycle boxplot:\n

       title        : plot title\n
       monthly_df   : pandas timeseries

    """
    # import required modules
    import numpy as np
    import matplotlib.pyplot as plt

    # new figure
    fig = plt.figure(figsize=[10, 8])
    ax = fig.add_subplot(111)

    # set axis
    if var == 1:
        ax.set_ylabel('Precipitación (mm)',
                      fontsize=12, fontweight='bold')
    else:
        ax.set_ylabel('Caudal (' + '$\mathregular{m^3}$' + '/s)',
                      fontsize=12, fontweight='bold')

    # plt annual cycle
    data = annual_cycle_df.values.reshape(12)
    errors = annual_cycle_error_df.values.reshape(12)
    index = np.arange(12)

    plt.bar(index, data, yerr=errors)

    plt.xticks(index, ['Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun', 'Jul', 'Ago',
                       'Sep', 'Oct', 'Nov', 'Dic'])

#    ax.set_xticks(1)
#    ax.set_xticklabels(['Ene', 'Feb', 'Mar', 'Abr', 'May', 'Jun', 'Jul', 'Ago', 'Sep', 'Nov',
#                        'Dic'])

    plt.title(title, fontweight='bold')

    return fig


# %% plot annual cycle boxplot
def plt_annual_cycle_boxplot(monthly_df, title, var):
    """Plot annual cycle boxplot:\n

       title        : plot title\n
       monthly_df   : pandas timeseries

    """
    # import required modules
    import matplotlib.pyplot as plt

    # new figure
    fig = plt.figure(figsize=[10, 8])
    ax = fig.add_subplot(111)

    # set axis
    ax.set_xlabel('Fecha', fontsize=12, fontweight='bold')
    if var == 1:
        ax.set_ylabel('Precipitación (mm)',
                      fontsize=12, fontweight='bold')
    else:
        ax.set_ylabel('Caudal (' + '$\mathregular{m^3}$' + '/s)',
                      fontsize=12, fontweight='bold')

    # plt annual cycle
    cols_list = list(monthly_df)[1:]
    monthly_df.boxplot(cols_list)
    plt.title(title, fontweight='bold')

    return fig


# %% plot precipitation mass curves
def plt_mass_curve(timeseries, title, var):
    """Plot pandas timeseries:\n

       title        : plot title\n
       timeseries   : pandas timeseries

    """
    # import required modules
    import matplotlib.pyplot as plt

    # new figure
    fig = plt.figure(figsize=[10, 6])
    ax = fig.add_subplot(111)

    # set axis
    ax.set_xlabel('Fecha', fontsize=12, fontweight='bold')
    if var == 1:
        ax.set_ylabel('Precipitación (mm)',
                      fontsize=12, fontweight='bold')
    else:
        ax.set_ylabel('Caudal (' + '$\mathregular{m^3}$' + '/s)',
                      fontsize=12, fontweight='bold')

    # plt timeseries
    timeseries = timeseries.cumsum()
    timeseries.plot()
    plt.title(title, fontweight='bold')

    return fig


def plt_hand_synthetic_rating_curve(depth, discharge, adjusted_depth, regression_parameters, fig_title):
    """
    Plot HAND synthetic rating curve.

    Parameters
    ----------
    depth: list
        Vector containing stream average depths
    discharge: list
        Vector containing stream hydraulic capacity
    adjusted_depth: list
        Adjusted depth from regression model
    regression_parameters: list
        Lineal regression parameters [slope, intercept, determination coefficient]
    fig_title: str
        Figure title

    Returns
    -------
    Figure containing HAND synthetic rating curve

    """
    # main imports
    import matplotlib.pyplot as plt

    # new figure
    fig, ax = plt.subplots(figsize=[10, 6])
    # plt.style.use(['classic'])

    # set axis
    ax.set_xlabel(r'$\mathbf{Streamflow\/[m^3/s]}$', fontsize=12, fontweight='bold')
    ax.set_ylabel('Stream depth [m]', fontsize=12, fontweight='bold')

    # plt scatter-plot and adjustment
    equation = r'$\mathrm{h = %0.3f \cdot Q^{%0.3f}}$' % (regression_parameters[1], regression_parameters[0])
    equation += '\n'
    equation += r'$\mathrm{R^2=%0.3f}$' % (regression_parameters[2])
    ax.plot(discharge, depth, label='Datos', marker='o', alpha=0.6)
    ax.plot(discharge, adjusted_depth, color='black', label='Ajuste', alpha=0.6)

    # legend and annotations
    handles, labels = ax.get_legend_handles_labels()
    handles = [handles[0], handles[1]]
    labels = [labels[0], labels[1]]
    legend = ax.legend(handles, labels, loc='upper left')
    plt.draw()
    p = legend.get_window_extent()
    ax.annotate(equation, (p.p0[0], p.p1[1]), (p.p0[0], p.p1[1] - 90), xycoords='figure pixels', zorder=9)

    # plot title
    plt.title(fig_title, fontweight='bold')

    return fig
