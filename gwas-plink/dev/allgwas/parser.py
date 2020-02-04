"""
functions to parse input

Date: 03/11/16
Author: Anja Gumpinger
"""
import logging
import os.path

"""
Configuration files
"""


class Config:
    """
    class of configuration files
    """

    def __init__(self, config_handle):
        logging.info('Parsing configuration file.')
        # every configuration file has input and output
        self.input = self._get_field(config_handle, 'data', 'input', None)
        self.output = self._get_field(config_handle, 'data', 'output', None)
        self.covar_fn = self._get_field(config_handle, 'data', 'covar_fn', None)
        self.pheno_fn = self._get_field(config_handle, 'data', 'pheno_fn', None)

    def _get_field(self, config_handle, level_1, level_2, default):
        """
        function reading from config file
        :param config_handle: handle to config file
        :param level_1: first level
        :param level_2: second level
        :param default: default value
        :return:
        """
        error_flag = False
        try:
            field = config_handle[level_1][level_2]
        except KeyError:
            # in case field does not exist in config file
            logging.info('%s, %s not specified, setting default: %s.' % (level_1, level_2, default))
            field = default
            error_flag = True

        if not field and not error_flag:
            logging.info('%s, %s not specified, setting default: %s.' % (level_1, level_2, default))
            field = default
        return field


class ConfigPlink(Config):
    """
    class of association configuration files, subclass of Config
    """

    def __init__(self, config_handle):
        Config.__init__(self, config_handle)

        # get plink parameters
        self.model = self._get_field(config_handle, 'assoc', 'model', 'assoc')
        self.covar_spec = self._get_field(config_handle, 'assoc', 'covar_spec', None)


class ConfigFastlmm(Config):
    """
    class of association configuration files, subclass of Config
    """

    def __init__(self, config_handle):
        Config.__init__(self, config_handle)

        # get fastlmm parameters
        self.k0 = self._get_field(config_handle, 'assoc', 'K0', None)
        self.k1 = self._get_field(config_handle, 'assoc', 'K1', None)


class ConfigFiltering(Config):
    """
    class of filtering configuration files, subclass of Config
    """

    def __init__(self, config_handle):
        Config.__init__(self, config_handle)
        self.hwe = self._get_field(config_handle, 'snps', 'hwe', float(1e-5))
        self.maf = self._get_field(config_handle, 'snps', 'maf', float(1e-2))


class ConfigGencon(Config):
    """
    class of genomic control configuration files, subclass of Config
    """

    def __init__(self, config_handle):
        Config.__init__(self, config_handle)
        self.model = self._get_field(config_handle, 'assoc', 'model', 'logistic')
        self.covar_spec = self._get_field(config_handle, 'assoc', 'covar_spec', None)
        self.num_pcs = self._get_field(config_handle, 'assoc', 'num_pcs', '10')


class ConfigEpistasis(Config):
    """
    class of genomic control configuration files, subclass of Config
    """

    def __init__(self, config_handle):
        Config.__init__(self, config_handle)
        self.epi1 = self._get_field(config_handle, 'epi', 'epi1', '0.0001')
        self.epi2 = self._get_field(config_handle, 'epi', 'epi2', '0.01')


class ConfigSetPlink(Config):
    """
    class of set-test configuration files for PLINK
    """

    def __init__(self, config_handle):
        Config.__init__(self, config_handle)
        self.model = self._get_field(config_handle, 'set', 'model', 'assoc')
        self.gene_file = self._get_field(config_handle, 'set', 'gene_file', None)
        self.mperm = self._get_field(config_handle, 'set', 'mperm', '1000')
        self.setp = self._get_field(config_handle, 'set', 'setp', '0.05')
        self.setr2 = self._get_field(config_handle, 'set', 'setr2', '0.5')
        self.setmax = self._get_field(config_handle, 'set', 'setmax', '5')
        self.margin = self._get_field(config_handle, 'set', 'margin', '20')


class ConfigSetFastlmm(Config):
    """
    """

    def __init__(self, config_handle):
        Config.__init__(self, config_handle)
        self.G0 = self._get_field(config_handle, 'set', 'G0', None)
        self.gene_file = self._get_field(config_handle, 'set', 'gene_file', None)
        self.test = self._get_field(config_handle, 'set', 'test', 'sc_davies')
        self.margin = self._get_field(config_handle, 'set', 'margin', '20')


"""
Text files
"""


class DelimFiles:
    """
    class of white-space delimited files
    """

    def __init__(self, filename=None):
        if filename:
            assert (os.path.isfile(filename)), 'File  %s does not exist.' % filename
            logging.info('File %s exists, reading.' % filename)
            self.__filename = filename
        else:
            logging.info('No filename specified.')

    @property
    def get_filename(self):
        return self.__filename

    def _read_col(self, col_idx, header=True):
        """
        read specified column
        :return:
        """
        lst = []
        with open(self.__filename, "r") as fin:
            if header:
                next(fin)
            for line in fin:
                parts = line.strip().split()
                lst.append(parts[col_idx])
        return lst

    def _convert(self, list, ctype):
        """
        convert all entries in a list to a specific type
        :param list: list containing elements to be converted
        :param ctype: type (float, int,...) elements should be converted to
        :return:
        """

        for idx, element in enumerate(list):
            try:
                list[idx] = ctype(list[idx])
            except ValueError:
                list[idx] = float('nan')
        return list


class PlinkLog:
    """
    class of PLINK log files
    """

    def __init__(self, filename):
        assert (os.path.isfile(filename)), 'File  %s does not exist.' % filename
        logging.info('File %s exists, reading.' % filename)
        self.__filename = filename
        self.__lambdagc = float(self._get_value(string="Genomic", split_idx=10)[:-1])

    @property
    def get_lambdagc(self):
        return self.__lambdagc

    def _get_value(self, string, split_idx):
        """
        extract one value from text file
        :param string: first word in line of interest
        :param split_idx: index of value of interest after splitting
        :return:
        """
        with open(self.__filename, "r") as fin:
            for line in fin:
                parts = line.strip().split()
                if not parts == []:
                    if string in parts:
                        return parts[split_idx]
                    else:
                        continue
        return None


"""
Specific text files
"""


class BimFiles(DelimFiles):
    """
    class of .bim files
    """

    def __init__(self, filename=None):
        DelimFiles.__init__(self, filename)
        self.__chrom = self._convert(self._read_col(col_idx=0, header=False), int)
        self.__snp = self._read_col(col_idx=1, header=False)
        self.__pos = self._convert(self._read_col(col_idx=3, header=False), int)

    @property
    def get_chromosomes(self):
        return self.__chrom

    @property
    def get_snps(self):
        return self.__snp

    @property
    def get_positions(self):
        return self.__pos


class FamFiles(DelimFiles):
    """
    class of .bim files
    """

    def __init__(self, filename=None):
        DelimFiles.__init__(self, filename)
        self.__fid = self._convert(self._read_col(col_idx=0, header=False), int)
        self.__iid = self._convert(self._read_col(col_idx=1, header=False), int)


    @property
    def get_family_ids(self):
        return self.__fid

    @property
    def get_individual_ids(self):
        return self.__iid


"""
Output Files
"""


class FlmmFiles(DelimFiles):
    """
    output of FaSTLMM analysis
    """

    def __init__(self, filename):
        DelimFiles.__init__(self, filename)
        self.__pvals = self._convert(self._read_col(col_idx=5), float)
        self.__chrom = self._convert(self._convert(self._read_col(col_idx=2), float), int)
        self.__pos = self._convert(self._read_col(col_idx=4), float)

    @property
    def get_pvalues(self):
        return self.__pvals

    @property
    def get_chromosome(self):
        return self.__chrom

    @property
    def get_position(self):
        return self.__pos


class PlinkFiles(DelimFiles):
    """
    class of standard plink-output file
    """

    def __init__(self, filename=None):
        if filename:
            DelimFiles.__init__(self, filename)
            self.__pvals = self._convert(self._read_col(col_idx=8), float)
            self.__stat = self._convert(self._read_col(col_idx=7), float)
            self.__chrom = self._convert(self._read_col(col_idx=0), int)
            self.__pos = self._convert(self._read_col(col_idx=2), int)

    @property
    def get_pvalues(self):
        return self.__pvals

    @property
    def get_statistics(self):
        return self.__stat

    @property
    def get_chromosome(self):
        return self.__chrom

    @property
    def get_position(self):
        return self.__pos


class PlinkAdjustedFiles(DelimFiles):
    """
    class of .adjusted PLINK files
    returns the chromosome, the snp-name, the raw and the gc-adjusted p-values
    """

    def __init__(self, filename=None):
        if filename:
            DelimFiles.__init__(self, filename)
            self.__gc = self._convert(self._read_col(col_idx=3), float)
            self.__unadj = self._convert(self._read_col(col_idx=2), float)
            self.__chrom = self._convert(self._read_col(col_idx=0), int)
            self.__snp = self._read_col(col_idx=1)

    @property
    def get_gc_pvals(self):
        return self.__gc

    @property
    def get_unadj_pvals(self):
        return self.__unadj

    @property
    def get_chromosome(self):
        return self.__chrom

    @property
    def get_snp(self):
        return self.__snp

    # def __getitem__(self, item):
    #     if isinstance(item, int):
    #         return self.__pvals[item]
    #     if isinstance(item, str):
    #         print "no"
    #         if isinstance(item, slice):
    #             return self.__pvals[slice]


class PlinkSetFiles(DelimFiles):
    """
    class of .set.mperm PLINK files
    returns the gene and the corresponding p-value, the snp-name, the raw and the gc-adjusted p-values
    """

    def __init__(self, filename=None):
        if filename:
            DelimFiles.__init__(self, filename)
            self.__pvalues = self._convert(self._read_col(col_idx=4), float)
            self.__genes = self._read_col(col_idx=0)

    @property
    def get_genes(self):
        return self.__genes

    @property
    def get_pvalues(self):
        return self.__pvalues
