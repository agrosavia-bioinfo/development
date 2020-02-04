"""
functions to execute plink with predefined list of command line parameters

Date: 02/11/16
Author: Anja Gumpinger
"""
import subprocess
import logging


def run(exec_fn, cmd_params):
    """
    run plink with command line parameters in cmd_params
    :param exec_fn: full path plink executable
    :type exec_fn: string

    :param cmd_params: command line parameters
    :type: list of strings

    :return:
    """

    p = subprocess.Popen([exec_fn] + cmd_params)
    return p
