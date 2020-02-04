"""
functions to do filtering of input

Date: 02/11/16
Author: Anja Gumpinger
"""
import allgwas.util
import logging


def filtering(configurations):
    """
    creates plink command line arguments to filter MAF and HWE
    :param configurations: object containing configurations
    :return:
    """
    # create log
    logging.info('Creating PLINK command line arguments to run filtering.')

    args = ["--file", configurations.input,
            "--make-bed",
            "--out", configurations.output,
            "--maf", configurations.maf,
            "--hwe", configurations.hwe]

    # create output directory, if it does not eist
    allgwas.util.check_dirs(configurations.output)

    return args


def eigendecomp(input_pref, output_pref, num_pcs):
    """
    creates plink command line arguments to compute eigendecomposition
    :param input_pref:
    :param output_pref:
    :param num_pcs:
    :return:
    """
    # create log
    logging.info('Creating PLINK command line arguments to run EVD.')

    args = ["--bfile", "%s" % input_pref,
            "--out", "%s" % output_pref,
            "--pca", "%s" % num_pcs]

    # create output directory, if it does not exist
    allgwas.util.check_dirs(output_pref)

    return args


def assoc(input_pref, output_pref, model, covar_fn=None, num_covar=None, pheno_fn=None):
    """
    creates command line arguments to run association testing
    :param input_pref: full path prefix to input file
    :param output_pref: full path prefix to output file
    :param model: which model to use (e.g. assoc, logistic, fisher, ...)
    :param covar_fn: full path filename of file containing covaritaes
    :param num_covar: number of covariates to use, format ("1-x" or "x, y, z" or mixture). if not specified, all
                      covariates in file will be used
    :return:
    """
    # create log
    logging.info('Creating PLINK command line arguments to run association analysis.')

    # generate vector containing command line arguments
    args = ["--noweb",
            "--bfile", "%s" % input_pref,
            "--out", "%s" % output_pref,
            "--%s" % model,
            "--adjust"]

    # if covariate file is specified, include in args
    if covar_fn:
        args = args + ["--covar", "%s" % covar_fn,
                       "--hide-covar"]
        if num_covar:
            args = args + ["--covar-number", num_covar]

    # if alternative phenotype file is specified
    if pheno_fn:
        args = args + ["--pheno", "%s" % pheno_fn,
                       "--allow-no-sex"]

    # create output directory, if it does not exist
    allgwas.util.check_dirs(output_pref)

    return args


def fast_epistasis(input_pref, output_pref, epi1, epi2):
    """
    creates command line arguments to run fast-epistasis

    :param input_pref: full path prefix to input files
    :param output_pref: full path prefix to output files
    :param epi1: maximum p-value for inclusion in main report
    :param epi2: maximum p-value to be counted as significant
    :return:
    """

    args = ["--bfile", "%s" % input_pref,
            "--out", "%s" % output_pref,
            "--fast-epistasis",
            "--epi1", "%s" % epi1,
            "--epi2", "%s" % epi2]

    return args


def epistasis(input_pref, output_pref, epi1, epi2):
    """
    creates command line arguments to run fast-epistasis

    :param input_pref: full path prefix to input files
    :param output_pref: full path prefix to output files
    :param epi1: maximum p-value for inclusion in main report
    :param epi2: maximum p-value to be counted as significant
    :return:
    """

    args = ["--bfile", "%s" % input_pref,
            "--out", "%s" % output_pref,
            "--epistasis",
            "--epi1", "%s" % epi1,
            "--epi2", "%s" % epi2]

    # create output directory, if it does not exist
    allgwas.util.check_dirs(output_pref)

    return args


def set_test(input_pref, output_pref, gene_list, model, mperm, setp, setmax, margin, setr2, pheno_fn=None):
    """
    creates command line arguments to run a set test

    :param input_pref: full path prefix to binary PLINK file
    :type input_pref: string

    :param output_pref: full path to output
    :type output_pref: string

    :param gene_list: full path of hg-file containing gene definitions
    :type gene_list: string

    :param model: model to use ('lrt' or 'sc_davies')
    :type model: string

    :param mperm: number of permutations
    :type mperm: string, int

    :param setp: minimum p-value for SNP to be included
    :type setp: string, float

    :param setmax: maximum number of SNPs per gene
    :type setmax: string, int

    :param margin: margin in [kb] around gene
    :type margin: string, int

    :param pheno_fn: full path to file containing phenotypes
    :type pheno_fn: string

    :param setr2:
    :type setr2:

    :param pheno_fn:
    :type pheno_fn:

    :return:
    """
    args = ["--bfile", "%s" % input_pref,
            "--out", "%s" % output_pref,
            "--%s" % model,
            "set-test",
            "--make-set", "%s" % gene_list,
            "--mperm", "%s" % mperm,
            "--set-p", "%s" % setp,
            "--set-max", "%s" % setmax,
            "--make-set-border", "%s" % margin,
            "--set-r2", "%s" % setr2]

    if pheno_fn:
        args = args + ["--pheno", "%s" % pheno_fn,
                       "--allow-no-sex"]

    # create output directory, if it does not exist
    allgwas.util.check_dirs(output_pref)

    return args