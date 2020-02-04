"""
function to generate the 2*4 figure plot containing manhattan and qq plot for different pop-struct correction methods

Date: 11/01/17
Author: Anja Gumpinger
"""
import scipy as sp
import allgwas.util


def preprocess_data(pvalue_lst, chromosome_lst, position_lst):
    """
    removes NaN-pvalues and orders SNPs in correct position according to pos_lst
    :param pvalue_lst:
    :param chromosome_lst:
    :param position_lst:
    :return:
    """

    # assert that all the input lists have the same length
    assert (len(chromosome_lst) == len(position_lst)), 'lists have to be of same length'
    assert (len(position_lst) == len(pvalue_lst)), 'lists have to be of same length'

    # remove SNPs that were assigned nan values
    nan_idx = sp.where(sp.isnan(sp.array(pvalue_lst)))[0]
    if len(nan_idx):
        # remove p-values with nan-values
        [pvalue_lst, chromosome_lst, position_lst] = allgwas.util.remove_list_entries(nan_idx,
                                                                                      [pvalue_lst, chromosome_lst,
                                                                                       position_lst])

        # make sure all nan entries were removed successfully
        assert (len(sp.where(sp.isnan(sp.array(pvalue_lst)))[0]) == 0), "Not all nan's were removed"

    # sort according to position
    chromosome_lst, position_lst, pvalue_lst = allgwas.util.sort_cpp_lists(chromosome_lst, position_lst, pvalue_lst)

    return chromosome_lst, position_lst, pvalue_lst


def make_plot(lst1, lst2, lst3, lst4, output_file):

    print 'test'