#!/usr/bin/env python3

# DFO-LS example: minimize the Rosenbrock function with noise
from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import dfols


def make_diagnostic_plot(df, n, xaxis_in_iters=False, show_restart_in_plots=True, obj_redn_in_tau=True, f0=None, fmin=None, max_ngradients=None, save_file=None):
    # df = soln.diagnostic_info (pandas DataFrame)
    # n = dimension of problem
    # xaxis_in_iters: True to plot against iteration of DFO-LS, False to plot against budget in evaluations
    # show_restart_in_plots: if using multiple restarts, then True plots vertical lines when each restart was called
    # obj_redn_in_tau: if True, plot objective reduction normalised so f0=1 and fmin=0 (need f0=starting objective value and fmin=true minimum objective value to be not None)
    # f0, fmin: see obj_redn_in_tau
    # max_ngradients: Maximum budget for DFO-LS (in units of simplex gradients, i.e. total evals/(n+1))
    # save_file: output image file; just display plot on screen if None
    
    font_size = 'small'

    ## These lines define values needed to do the plot of iteration types ##
    # from diagnostic_info.py
    ITER_VERY_SUCCESSFUL = "Very successful"  # ratio >= 0.7, no geometry update
    ITER_SUCCESSFUL = "Successful"  # 0.1 <= ratio < 0.7, no geometry update
    ITER_ACCEPTABLE_GEOM = "Acceptable (geom fixed)"  # 0 <= ratio < 0.1, with geometry update
    ITER_ACCEPTABLE_NO_GEOM = "Acceptable (geom not fixed)"  # 0 <= ratio < 0.1, without geometry update
    ITER_UNSUCCESSFUL_GEOM = "Unsuccessful (geom fixed)"  # ratio < 0, with geometry update
    ITER_UNSUCCESSFUL_NO_GEOM = "Unsuccessful (geom not fixed)"  # ratio < 0, without geometry update (possibly rho reduced)
    ITER_SAFETY = "Safety"  # safety step taken (||s|| too small compared to rho)

    # Info for plotting
    ITER_PLOTTING_INFO = {}
    ITER_PLOTTING_INFO[ITER_VERY_SUCCESSFUL] = 5  # value to plot, text
    ITER_PLOTTING_INFO[ITER_SUCCESSFUL] = 4  # value to plot, text
    ITER_PLOTTING_INFO[ITER_ACCEPTABLE_NO_GEOM] = 4  # value to plot, text
    ITER_PLOTTING_INFO[ITER_ACCEPTABLE_GEOM] = 3  # value to plot, text
    ITER_PLOTTING_INFO[ITER_UNSUCCESSFUL_GEOM] = 2  # value to plot, text
    ITER_PLOTTING_INFO[ITER_UNSUCCESSFUL_NO_GEOM] = 1  # value to plot, text
    ITER_PLOTTING_INFO[ITER_SAFETY] = 0  # value to plot, text
    ITER_PLOTTING_INFO["n/a"] = -1
    ITER_PLOTTING_INFO[None] = -1

    ITER_PLOT_LABELS = {-1: "", 0: "Safety", 1: "Bad", 2: "Bad + G", 3: "Ok + G", 4: "Ok", 5: "Good"}
    ## End iteration type info ##

    # Where to plot vertical lines (from restarts)
    nruns = df['nruns'].as_matrix()
    skips = nruns[1:] != nruns[:-1]
    idx_skips = np.where(skips)[0]
    budget_where_restarts_occurred = list(df['nf'][idx_skips].as_matrix() / float(n + 1))
    
    vlines = []
    if xaxis_in_iters:
        for i in range(len(idx_skips)):
            vlines.append(idx_skips[i] + 0.5)
    else:
        for b in budget_where_restarts_occurred:
            vlines.append(b - 0.5 / (n + 1))

    plt.ioff()  # non-interactive mode
    plt.figure(1, figsize=(9, 6))  # size in inches (usually save 100dpi, so pixels is multiplied by 100)

    xaxis = (df['iters_total'] if xaxis_in_iters else df['nf'] / (n + 1)).as_matrix()
    xmin = np.min(xaxis)
    xmax = np.max(xaxis)
    if not xaxis_in_iters and max_ngradients is not None:
        xmax = min(xmax, max_ngradients)  # normalise x-axis by total budget for easy comparison
    xlabel = "Iteration" if xaxis_in_iters else "Budget (gradients)"

    plt.subplot(2, 2, 1)  # (2x2) grid, first subplot
    ax = plt.gca()  # current axes
    ax.semilogy(xaxis, df['delta'], color='b', linestyle='-', linewidth=2.0, label=r"$\Delta_k$")
    ax.semilogy(xaxis, df['rho'], color='b', linestyle='--', linewidth=1, label=r"$\rho_k$")
    ax.semilogy(xaxis, df['norm_gk'], color='r', linestyle='-', linewidth=2.0, label=r"$||g_k||$")
    ax.semilogy(xaxis, df['norm_sk'], color='g', linestyle='-', linewidth=2.0, label=r"$||s_k||$")
    if 'poisedness' in df:
        ax.semilogy(xaxis, df['poisedness'], color='c', linestyle='-', linewidth=2.0, label="Poisedness")
    if show_restart_in_plots:
        for v in vlines:
            ax.axvline(x=v, linewidth=2.0, linestyle='-', color='k')
    leg = ax.legend(loc='upper right', fontsize=font_size, fancybox=True)
    ax.set_xlim(xmin, xmax)
    # ax.set_xlabel(xlabel)
    ax.set_title('Convergence info')
    ax.grid()

    plt.subplot(2, 2, 2)  # (2x2) grid, second subplot
    ax = plt.gca()  # current axes
    iter_types = [ITER_PLOTTING_INFO[t] for t in df['iter_type'].tolist()]
    ax.plot(list(xaxis), iter_types, color='b', linestyle='none', linewidth=2.0, marker='.',
            markersize=10, label="Iteration type")
    if show_restart_in_plots:
        for v in vlines:
            ax.axvline(x=v, linewidth=2.0, linestyle='-', color='k')
    tick_levels = []
    tick_labels = []
    for t in ITER_PLOTTING_INFO:
        if t != "n/a" and t is not None:  # and t != ITER_INIT:
            tick_levels.append(ITER_PLOTTING_INFO[t])
            tick_labels.append(ITER_PLOT_LABELS[ITER_PLOTTING_INFO[t]])
    # ax.axis([df['iters_total'].min(), df['iters_total'].max(), min(tick_levels) - 0.1, max(tick_levels) + 0.1])  # (xlow, xhigh, ylow, yhigh)
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(min(tick_levels) - 0.1, max(tick_levels) + 0.1)
    ax.set_yticks(tick_levels)
    ax.set_yticklabels(tick_labels)
    # ax.set_xlabel(xlabel)
    ax.set_title("Iteration type")
    ax.grid()

    plt.subplot(2, 2, 3)  # (2x2) grid, third subplot
    ax = plt.gca()  # current axes
    ax.semilogy(xaxis, df['interpolation_change_J_norm'], color='b', linestyle='-', linewidth=2.0, label=r"$||J_k-J_{k-1}||_2$")
    if show_restart_in_plots:
        for v in vlines:
            ax.axvline(x=v, linewidth=2.0, linestyle='-', color='k')
    ax.set_xlabel(xlabel)
    ax.set_xlim(xmin, xmax)
    leg = ax.legend(loc='lower right', fontsize=font_size, fancybox=True)
    ax.set_title('Jacobian info')
    ax.grid()
    
    plt.subplot(2, 2, 4)  # (2x2) grid, fourth subplot
    ax = plt.gca()  # current axes
    fopts = np.maximum(df['fk'].as_matrix(), 1e-12)  # floor values at 1e-12
    if obj_redn_in_tau:
        fopts = (fopts - fmin) / (f0 - fmin)  # switch to tau format
    ax.semilogy(xaxis, fopts, color='b', linestyle='-', linewidth=2.0, label="Objective reduction")
    if not obj_redn_in_tau:
        ax.axhline(y=max(fmin, 1e-12), linewidth=2.0, linestyle='--', color='r')
    if show_restart_in_plots:
        for v in vlines:
            ax.axvline(x=v, linewidth=2.0, linestyle='-', color='k')
    ax.set_xlabel(xlabel)
    ax.set_xlim(xmin, xmax)
    # leg = ax.legend(loc='lower right', fontsize=font_size, fancybox=True)
    ax.set_title('Normalised obj redn' if obj_redn_in_tau else 'Objective reduction')
    ax.grid()

    if save_file is not None:
        plt.savefig(save_file, bbox_inches='tight')
    else:
        plt.show()
    return


if __name__ == '__main__':
    # Define the objective function
    def rosenbrock(x):
        return np.array([10.0 * (x[1] - x[0] ** 2), 1.0 - x[0]])
    
    # Modified objective function with additive Gaussian noise
    def rosenbrock_noisy(x, sigma=1e-2):
        return rosenbrock(x) + sigma * np.random.normal(size=(2,))
    
    # Define the starting point
    x0 = np.array([-1.2, 1.0])
    
    # Set random seed (for reproducibility)
    np.random.seed(0)
    
    # Extra parameters to save diagnostic info
    params = {}
    params["logging.save_diagnostic_info"] = True
    params["logging.save_xk"] = False  # optional: save best x value at every iteration
    params["logging.save_rk"] = False  # optional: save r(x) corresponding to best x value at every iteration (automatically save overall cost function f(xk))
    params["logging.save_poisedness"] = True  # optional: save 'poisedness' (measure of geometry of interpolation set: smaller is better); this can be expensive to calculate if there are many parameters to be optimised over
    
    # Call DFO-LS
    soln = dfols.solve(rosenbrock_noisy, x0, objfun_has_noise=True, user_params=params)
    
    # Display output
    print(soln)
    
    # Get diagnosic info
    df = soln.diagnostic_info  # pandas DataFrame object, one row per iteration (see Github documentation for what each column represents)
    # save to df csv using: df.to_csv('dfols_info.csv', index=False)
    # Read from csv: import pandas as pd; df = pd.read_csv('dfols_info.csv')
    
    # True f0, fmin for Rosenbrock problem
    n = len(x0)
    r0 = rosenbrock(x0)
    f0 = np.dot(r0, r0)
    fmin = 0.0
    
    save_file = None  # display plot on screen
    #save_file = 'dfols_results.png'  # save plot to file, can change filetype if desired
    make_diagnostic_plot(df, n, xaxis_in_iters=False, show_restart_in_plots=True, obj_redn_in_tau=False, f0=f0, fmin=fmin, save_file=save_file)
    
    print("Done")
