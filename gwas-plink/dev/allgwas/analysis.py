"""
all functions to analyse GWAS results

Date: 04/11/16
Author: Anja Gumpinger
"""
import allgwas.parser
import scipy as sp
import matplotlib.pyplot as plt
import math
import allgwas.util
import logging


def qq(obs_pval, lambda_gc=None, output_fn=None, show_labels=True, x_limits=None, y_limits=None):
    """
    QQ-plot of p-values

    :param obs_pval: observed pvalues
    :type obs_pval: list

    :param lambda_gc: genomic inflation
    :type: float

    :param output_fn: full path of file to store figure in
    :type output_fn: string

    :param show_labels: whether or not to display axis labels
    :type show_labels: boolean

    :param x_limits: range (minimum and maximum) of x-axis
    :type x_limits: list of float

    :param y_limits: range (minimum and maximum) of y-axis
    :type y_limits: list of float

    :return:
    """
    # log info
    logging.info('Generating QQ-plot.')

    # generate expected p-values
    exp_pval = sp.linspace(0, 1, len(obs_pval) + 1)[1:].tolist()

    # check if p-values are nans, remove them
    nan_idx = sp.where(sp.isnan(sp.array(obs_pval)))[0]
    if len(nan_idx):
        logging.warning('Of %s analyzed SNPs, %s have \'nan\' p-values. Removing Snps.'
                        % (len(obs_pval), len(nan_idx)))
        # remove p-values with nan-values
        [obs_pval, exp_pval] = allgwas.util.remove_list_entries(nan_idx, [obs_pval, exp_pval])

        # make sure all nan entries were removed successfully
        assert (len(sp.where(sp.isnan(sp.array(obs_pval)))[0]) == 0), "Not all nan's were removed"
        logging.info('Nans removed.')

    # transform to -log10
    obs_pval = [-sp.log10(p) for p in obs_pval]
    exp_pval = [-sp.log10(p) for p in exp_pval]

    # sort
    obs_pval_sorted = sorted(obs_pval)
    exp_pval_sorted = sorted(exp_pval)

    # get maximum value of p-values (to get limits of plot properly)
    max_val = math.ceil(max([obs_pval_sorted[-1], exp_pval[-2]]) + 1)

    # generate QQ-plot (observed vs. expected p-values, and line indicating uniform distribution)
    fig, ax = plt.subplots()
    ax.plot([0, max_val], [0, max_val], '-k')
    ax.plot(exp_pval_sorted, obs_pval_sorted, '.', color="#1F78B4")

    if lambda_gc:
        plt.title("inflation: %.3s" % lambda_gc)

    # add grid and axis labels
    plt.grid(True)
    if show_labels:
        ax.set_xlabel('expected -log10(pvalues)')
        ax.set_ylabel('observed -log10(pvalues)')

    # change the tick size
    ax.xaxis.set_tick_params(labelsize=15)
    ax.yaxis.set_tick_params(labelsize=15)

    # if specified: set axis limits:
    if x_limits:
        ax.set_xlim(x_limits)
    if y_limits:
        ax.set_ylim(y_limits)

    # save or show figure, depending on specification
    if output_fn:
        print "saving q-q plot to %s" % output_fn
        # check if output folder exists and save
        allgwas.util.check_dirs(filename=output_fn)
        plt.savefig(output_fn, bbox_inches='tight')
    else:
        plt.show()
    pass


def manhattan(chromosome_lst, position_lst, pvalue_lst, output_fn=None, show_labels=True, y_limits=None):
    """
    Manhattan plot of p-values

    :param chromosome_lst: list containing chromosomes
    :type chromosome_lst: list of int

    :param position_lst: list containing SNP positions
    :type position_lst: list of int

    :param pvalue_lst: list containing p-values
    :type pvalue_lst: list of float

    :param output_fn: filename of figure to generate
    :type output_fn: string

    :param show_labels: whether or not to display axis labels
    :type show_labels: boolean

    :param y_limits: range (minimum and maximum) of y-axis
    :type y_limits: list of float
    :return:
    """

    # assert that all the input lists have the same length
    assert (len(chromosome_lst) == len(position_lst)), 'lists have to be of same length'
    assert (len(position_lst) == len(pvalue_lst)), 'lists have to be of same length'

    # remove SNPs that were assigned nan values
    nan_idx = sp.where(sp.isnan(sp.array(pvalue_lst)))[0]
    if len(nan_idx):
        logging.warning('Of %s analyzed SNPs, %s have \'nan\' p-values. Removing Snps.'
                        % (len(pvalue_lst), len(nan_idx)))
        # remove p-values with nan-values
        [pvalue_lst, chromosome_lst, position_lst] = allgwas.util.remove_list_entries(nan_idx,
                                                                                      [pvalue_lst, chromosome_lst,
                                                                                       position_lst])

        # make sure all nan entries were removed successfully
        assert (len(sp.where(sp.isnan(sp.array(pvalue_lst)))[0]) == 0), "Not all nan's were removed"

    # sort lists by chromosome and position
    chromosome_lst, position_lst, pvalue_lst = allgwas.util.sort_cpp_lists(chromosome_lst, position_lst, pvalue_lst)

    # transform p-values to -log10
    pv_arr = sp.asarray([-sp.log10(p) for p in pvalue_lst])

    # generate color vector (each entry corresponds to one chromosome)
    col_vec = ['#A6CEE3' if i % 2 == 0 else '#1F78B4' for i in sp.unique(chromosome_lst)]

    # generate manhattan plot
    fig, ax = plt.subplots()

    # initialize location and labels of x-axis ticks
    xtick_loc = []
    xtick_label = []

    # compute Bonferroni correction
    bonf = 7.2e-8
    log_bonf = -sp.log10(bonf)

    # plot p-values versus SNPs
    for chrom in sp.unique(chromosome_lst):
        idx = sp.asarray(sp.where(chromosome_lst == chrom)[0])
        ax.plot([idx[0], idx[-1]], [log_bonf, log_bonf], '-r', lw=1.5)
        ax.plot(idx, pv_arr[idx], '.', color=col_vec[chrom - 1], markersize=4)

        # generate locations for the x-axis ticks and resp. labels
        xtick_loc.append(idx[0] + (idx[-1] - idx[0]) / 2)
        xtick_label.append('%s' % chrom)

    # set x-axis limits
    ax.set_xlim([0, idx[-1]])

    # if specified, set y-limits:
    if y_limits:
        print y_limits
        ax.set_ylim(y_limits)

    # add grid, add labels for x-axis
    ax.get_xaxis().set_ticks(xtick_loc)
    ax.get_xaxis().set_ticklabels(xtick_label, fontsize=15, rotation=90)
    plt.gca().xaxis.grid(False)
    plt.gca().yaxis.grid(True)

    # increase the fontsize of the y-ticks
    ax.yaxis.set_tick_params(labelsize=15)

    # axis labels
    if show_labels:
        ax.set_xlabel('chromosome')
        ax.set_ylabel('significance [-log10(pvalue)]')

    if output_fn:
        # check if output folder exists and save
        allgwas.util.check_dirs(filename=output_fn)
        plt.savefig(output_fn, bbox_inches='tight')
    else:
        plt.show()
    pass


def scree_plot(lambda_gc_vec, figure_fn=None):
    """
    Scree plot of lambda_gc values

    :param lambda_gc_vec: list containing lambda_gc values
    :type lambda_gc_vec: list of float

    :param figure_fn: filename of output figure
    :type figure_fn: string

    :return:
    """

    # log info
    logging.info('Generating scree plot.')

    # generate scree plot
    fig, ax = plt.subplots()

    # set locations on x-axis
    x_pos = sp.arange(0, len(lambda_gc_vec))

    # dot-line plot of genomic inflation
    ax.plot(x_pos, lambda_gc_vec, '-', color='#1F78B4')
    ax.plot(x_pos, lambda_gc_vec, '.', color='#1F78B4', markersize=8)

    # plot star for lambda_gc with minimal distance to 1
    md_idx = min(range(len(lambda_gc_vec)), key=lambda i: abs(lambda_gc_vec[i]-1))
    ax.plot(x_pos[md_idx], lambda_gc_vec[md_idx], '*', color='r', markersize=12, label='best')

    # set axis labels
    ax.set_xlabel('number of leading PCs used as covariates')
    ax.set_ylabel('genomic inflation factor')

    # set grid
    plt.grid(True)
    plt.legend(numpoints=1)

    # change the tick size
    ax.xaxis.set_tick_params(labelsize=15)
    ax.yaxis.set_tick_params(labelsize=15)

    # save
    if figure_fn:
        # check if output folder exists and save
        allgwas.util.check_dirs(filename=figure_fn)
        plt.savefig(figure_fn, bbox_inches='tight')
    else:
        plt.show()
    pass

