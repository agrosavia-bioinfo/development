"""
All functions to convert files from one format to another

Date: 09/01/17
Author: Anja Gumpinger
"""
import parser


def parse_tair(input_file, output_file):
    """
    function to parse TAIR gff file

    :param input_file: full path filename of TAIR file
    :type input_file: string

    :param output_file: full path filename of parsed file
    :type output_file: string

    :return:
    """

    # open output file to write to
    fout = open(output_file, 'w')

    # open input file in read mode
    with open(input_file, 'r') as fin:
        for line in fin:

            # split the line in its single fragments
            parts = line.strip().split()

            # if the line contains gene information
            if parts[2] == 'gene':

                # extract the chromosome, the start and end position of the gene, and the gene ID
                chrom = parts[0].split('Chr')[1]
                start = parts[3]
                end = parts[4]
                id = parts[-1].split(';')[0].split('=')[1]

                try:
                    # check if chromosome identifier can be converted to int (to exclude chromosome C and M)
                    int(chrom)
                    # write this to output file
                    fout.write("%s\n" % '\t'.join([chrom, start, end, id]))
                except ValueError:
                    continue
    fout.close()

    pass


def parse_tair_interaction(input_file, output_file):
    """
    function to parse TAIR interaction file

    :param input_file: full path filename to TAIR interaction file
    :type input_file: string

    :param output_file: full path filename of parsed file
    :type output_file: string
    :return:
    """

    # open output file to write to
    fout = open(output_file, 'w')

    # many interactions occur multiple times; to skip them, always remember previous locus
    locus1_prev = ''
    locus2_prev = ''

    # open input file in read mode
    with open(input_file, 'r') as fin:
        next(fin)
        for line in fin:

            # split the line in its single fragments
            parts = line.strip().split()

            # extract the two interacting loci
            locus1 = parts[0]
            locus2 = parts[2]

            if not '\t'.join([locus1, locus2]) == '\t'.join([locus1_prev, locus2_prev]):
                fout.write('%s\n' % '\t'.join([locus1, locus2]))

            locus1_prev = locus1
            locus2_prev = locus2
    fout.close()


def parse_eigensoft_eigenvectors(eigensoft_file, fam_file, output_file):
    """
    problem with files from eigensoft: they only contain individual ID. This function matches the individual ID with
    its family ID from a fam-file and writes the eigenvectors in PLINK's covariate/phenotype file format

    :param eigensoft_file: {filename}.evec file generated with EIGENSOFT's smartpca.perl function
    :type eigensoft_file: string

    :param fam_file: full path to fam-file containing the individual information
    :type fam_file: string

    :param output_file: full path filename of output file
    :type output_file: string

    :return:
    """

    # parse input file
    fam_file = parser.FamFiles(fam_file)
    iid_vec = fam_file.get_individual_ids
    fid_vec = fam_file.get_family_ids

    # create dict to allow for faster searching through individuals
    id_dict = {}
    for (iid, fid) in zip(iid_vec, fid_vec):
        id_dict[iid] = fid

    # open and write to output file
    fout = open(output_file, 'w')

    with open(eigensoft_file, 'r') as fin:
        next(fin)
        for line in fin:
            parts = line.strip().split()
            fout.write("%s\t%s\n" % (id_dict[int(parts[0])], '\t'.join(parts[:-1])))
    fout.close()

    pass


def parse_gene_pvalues(pvalue_file, output_file):
    """
    extract pvalues from PLINK set-test output file, and bring it into two columned layout for dmGWAS

    :param pvalue_file: full path to file containing the p-values that were generated with PLINK's --set-test
    :type pvalue_file: string

    :param output_file: full path filename of output file
    :type output_file: string
    """

    # parse the input file
    plink_file = parser.PlinkSetFiles(pvalue_file)
    pval_vec = plink_file.get_pvalues
    gene_vec = plink_file.get_genes

    # open and write to output file
    fout = open(output_file, 'w')

    for (gene, pv) in zip(gene_vec, pval_vec):
        # dmGWAS cannot deal with p-values == 1/0, so replace them
        if pv == 1.0:
            pv = 0.9999
        if pv == 0.0:
            pv = 0.0001
        fout.write("%s\t%s\n" % (gene, pv))
    fout.close()

    pass
