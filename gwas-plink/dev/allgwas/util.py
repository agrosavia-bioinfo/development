"""
helping functions

Date: 04/11/16
Author: Anja Gumpinger
"""
import os
import allgwas.plink
import subprocess


def check_dirs(filename):
    """
    check if the path of a file exists, if not, create it

    :param filename: full path of file for which path should be checked
    :type filename: string

    :return:
    """
    if not os.path.isdir(os.path.dirname(filename)):
        os.makedirs(os.path.dirname(filename))

    pass


def remove_list_entries(indices, lists):
    """
    removes entries in 'indices' fro all lists in 'lists'

    :param indices: indices of elements to remove
    :type indices: list of indices

    :param lists: list of lists, or single list from which entries with indices 'indices' should be removed
    :type lists: list

    :return:
    """
    # check if dealing with list of lists, or with one single list
    if isinstance(lists[0], list):
        new_lists = []
        for lst in lists:
            new_lists.append([i for j, i in enumerate(lst) if j not in indices])

    # if single list
    else:
        new_lists = [i for j, i in enumerate(lists) if j not in indices]

    return new_lists


def create_pheno_file(fam_file_prefix):
    """
    extracts the phenotype from a fam-file (6th column) and writes it to a separate file

    :param fam_file_prefix: full path prefix of PLINK fam file
    :type fam_file_prefix: string

    :return:
    """
    # generate and open output file
    pheno_fn = "%s_pheno.txt" % fam_file_prefix
    fout = open(pheno_fn, "w")

    with open("%s.fam" % fam_file_prefix, 'r') as fin:
        for line in fin:

            # split line into its fragments
            parts = line.strip().split()

            # write fid, iid and phenotype to output file
            fout.write("%s\t%s\t%s\n" % (parts[0], parts[1], parts[5]))
    fout.close()
    return pheno_fn


def print_configurations(config):
    """
    prints settings in configuration

    :param config: contains all information required to run script
    :type config: instance of class config

    :return:
    """
    print "Settings used"
    for key in config.__dict__.keys():
        print "%s:\t %s" % (key, config.__dict__[key])

    pass


def recode_phenotype(input_fn, output_fn, original, recoded):
    """
    recode a phenotype; PLINK has the standart assignment:
    -9  missing
    0   missing
    1   unaffected
    2   affected

    :param input_fn: full path to input file (to be recoded)
    :type input_fn: string

    :param output_fn: full path to output file (recoded)
    :type output_fn: string

    :param original: original value
    :type original: int

    :param recoded: replacement value
    :type recoded: int
    :return:
    """
    # open output file
    fout = open(output_fn, 'w')

    with open(input_fn, "r") as fin:
        next(fin)
        for line in fin:

            # split line into its fragments
            parts = line.strip().split()

            # if value to be recoded is encountered
            if int(parts[2]) == original:
                
                # recode
                parts[2] = str(recoded)

            # write to output file
            fout.write("%s\n" % "\t".join(parts))
    pass


def create_fastlmm_set_file(snp_file, gene_annotation_file, margin, output_fn, pheno_fn=None):
    """
    from hg files (downloaded from PLINK homepage), generate gene-files for running fast-LMM

    :param snp_file: prefix to binary PLINK file containing SNPs to be mapped to genes
    :type snp_file: string

    :param gene_annotation_file: full path to file of hg-built of gene-positions
    :type gene_annotation_file: string

    :param margin: size [kb] of margin around genes. SNPs lying in margin will also be assigned to gene
    :type margin: string

    :param output_fn: full path filename of output
    :type output_fn: string

    :param pheno_fn: full path filename to phenotype-file (needed to run association, just as a dummy), catgorical!!!
    :type pheno_fn: string

    :return: full path of file containing SNP-gene mapping
    """
    # run plink set test to obtain mapping
    cmd_assoc = allgwas.plink.cmd.set_test(input_pref=snp_file,
                                           output_pref="%s_tmp" % snp_file,
                                           model='assoc',
                                           gene_list=gene_annotation_file,
                                           mperm=1,
                                           setp=1,
                                           setmax=1000,
                                           margin=margin,
                                           setr2=1,
                                           pheno_fn=pheno_fn)

    # surpress console output
    cmd_assoc = cmd_assoc + ["--silent"]

    # run plink
    proc = allgwas.plink.execute.run(exec_fn='plink', cmd_params=cmd_assoc)
    proc.wait()

    # extract mapping from PLINK file
    fout = open("%s" % output_fn, "w")

    # write header
    fout.write('snp\tset\n')

    with open("%s_tmp.assoc.set.mperm" % snp_file, "r") as fin:
        next(fin)
        for line in fin:
            parts = line.strip().split()

            if not parts[5] == 'NA':
                snps = parts[5].split('|')
                for snp in snps:
                    fout.write("%s\t%s\n" % (snp, parts[0]))
    fout.close()

    # delete intermediate files (not needed)
    subprocess.call(["rm", "%s_tmp.assoc" % snp_file])
    subprocess.call(["rm", "%s_tmp.assoc.set.mperm" % snp_file])
    subprocess.call(["rm", "%s_tmp.log" % snp_file])
    subprocess.call(["rm", "%s_tmp.nosex" % snp_file])

    pass


def create_plink_set_file(snp_file, gene_annotation_file, margin, output_fn, pheno_fn=None):
    """
    from hg files (downloaded from PLINK homepage), generate gene-files for running fast-LMM

    :param snp_file: prefix to binary PLINK file containing SNPs to be mapped to genes
    :type snp_file: string

    :param gene_annotation_file: full path to file of hg-built of gene-positions
    :type gene_annotation_file: string

    :param margin: size [kb] of margin around genes. SNPs lying in margin will also be assigned to gene
    :type margin: string

    :param output_fn: full path filename of output
    :type output_fn: string

    :param pheno_fn: full path filename to phenotype-file (needed to run association, just as a dummy), catgorical!!!
    :type pheno_fn: string

    :return: full path of file containing SNP-gene mapping
    """

    # run plink set test to obtain mapping
    cmd_assoc = allgwas.plink.cmd.set_test(input_pref=snp_file,
                                           output_pref="%s_tmp" % snp_file,
                                           model='assoc',
                                           gene_list=gene_annotation_file,
                                           mperm=1,
                                           setp=1,
                                           setmax=1000,
                                           margin=margin,
                                           setr2=1,
                                           pheno_fn=pheno_fn)

    # surpress console output
    cmd_assoc = cmd_assoc + ["--silent"]

    # run plink
    proc = allgwas.plink.execute.run(exec_fn='plink', cmd_params=cmd_assoc)
    proc.wait()

    # extract mapping from PLINK file
    fout = open("%s" % output_fn, "w")

    with open("%s_tmp.assoc.set.mperm" % snp_file, "r") as fin:
        next(fin)
        for line in fin:

            # extract fragments
            parts = line.strip().split()

            if not parts[5] == 'NA':

                # get SNPs
                snps = parts[5].split('|')

                # write set name
                fout.write("%s\n" % parts[0])

                # write all SNPs
                for snp in snps:
                    fout.write("%s\n" % snp)

                # write 'END'
                fout.write('END\n\n')

    fout.close()

    # delete intermediate files (not needed)
    subprocess.call(["rm", "%s_tmp.assoc" % snp_file])
    subprocess.call(["rm", "%s_tmp.assoc.set.mperm" % snp_file])
    subprocess.call(["rm", "%s_tmp.log" % snp_file])
    subprocess.call(["rm", "%s_tmp.nosex" % snp_file])

    pass


def sort_cpp_lists(chromosome_lst, position_lst, pvalue_lst):
    """
    sorts the chromosome, snp-position and p-value list, such that all are sorted by chromosomes, and each chromosome is
    sorted by the snp-position

    :param chromosome_lst: chromosomes
    :type chromosome_lst: list of int

    :param position_lst: SNP positions
    :type position_lst: list of int

    :param pvalue_lst: p-values
    :type pvalue_lst: list of float

    :return: sorted lists
    """

    # 0. Step: zip lists
    zipped = zip(chromosome_lst, position_lst, pvalue_lst)

    # 1. Step: sort them by position
    zipped.sort(key=lambda t: t[1])

    # 2. Step: sort by chromosome
    zipped.sort(key=lambda t: t[0])

    # 3. Step: unzip
    unzipped = zip(*zipped)

    # 4. Step: return lists (chromosome, position, pvalue)
    return list(unzipped[0]), list(unzipped[1]), list(unzipped[2])


def get_position_from_bim(snp_lst, bim_file):
    """
    returns the positions of every snp in snp_lst, read from bim_file

    :param snp_lst: list of SNPs for which position should be returned
    :type snp_lst: list

    :param bim_file: bim_file containing SNP information
    :type bim_file: string

    :return: list of SNP positions
    """

    # create instance of BimFiles and get snps and positions
    bim = allgwas.parser.BimFiles(bim_file)
    bim_snps = bim.get_snps
    bim_positions = bim.get_positions

    # read them into dict (easier to search)
    snp_dict = {}
    while bim_snps:
        snp_dict[bim_snps.pop()] = bim_positions.pop()

    # create list of positions to return
    positions = []
    for snp in snp_lst:
        positions.append(snp_dict[snp])

    return positions








