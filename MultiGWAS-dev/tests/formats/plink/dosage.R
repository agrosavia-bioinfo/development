# will use data.table for faster I/O, but it's not strictly necessary for this script to work
library(data.table)

# function to perform conversion
# call outside of R with plink2 (latest build, 2019-09 or later)
# infile is a file of N samples on rows and P SNPs on columns, a header, and leftmost column with sample IDs
# remaining contents of infile are hard-called dosages (0,1,2)
# call PLINK on result of files like this:
# ~/bin/plink2 --import-dosage $OUT_DOSAGE noheader --fam $OUT_FAM --make-bed --out $YOUR_OUTFILE_PREFIX

convert.raw.to.plink2 = function(infile, out.dosage, out.fam){

    # load genotype dosages
    geno = fread(infile)

    # transpose the allele dosages
    # use as.matrix to ensure that we (efficiently) transpose the allele dosage numbers and nothing else
    tgeno = t(as.matrix(geno))

    # build the dosage file itself
    geno.dosage = data.table(cbind(paste0("snp", 1:ncol(geno)), "A", "T", tgeno))

    # could add header, this one should work with --id-delim "-" in PLINK call
    #colnames(ceu.dosage) = c("SNP", "A1", "A2", paste0("id", 1:nrow(ceu), "-id", 1:nrow(ceu)))

    # write the dosage file to disk
    # set col.names=T if you want a header on the file  
    fwrite(geno.dosage, file = out.dosage, quote = F, sep = "\t", col.names=F)

    # build a dummy fam file and write that to disk
    geno.fam = data.table(cbind(paste0("id", 1:nrow(geno)), paste0("id", 1:nrow(geno)), 0, 0, 0, -9))
    fwrite(geno.fam, file = out.fam, quote = F, sep = "\t", col.names=F)

    return()

}

infile = "plink.raw"
out.dosage = "dplink.dosage"
out.fam     = "dplink.fam"

convert.raw.to.plink2 (infile, out.dosage, out.fam)
