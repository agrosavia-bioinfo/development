# Parameters
data_dir=/home/gwasuser/data/a.thaliana

# ----------------------------------------------------------------
# 4.1
# ----------------------------------------------------------------
#folder=4.1_data_formats

#script=example_4.1_binary2plain.sh
#echo "> $script"
#time $EXAMPLES/${folder}/${script} > /dev/null 2>&1

#script=example_4.1_plain2binary.sh
#echo "> $script"
#time $EXAMPLES/${folder}/${script} > /dev/null 2>&1
## Remove the files created by these scripts
#rm $data_dir/atwell2010_geno_new.*

# ----------------------------------------------------------------
# 4.2 
# ----------------------------------------------------------------
#folder=4.2_data_preprocessing

#script=example_4.2_preprocessing.sh
#echo "> $script"
#time $EXAMPLES/${folder}/${script} > /dev/null 2>&1

# ----------------------------------------------------------------
# 4.3 
# ----------------------------------------------------------------
out_dir=/home/gwasuser/output/univariate_gwas/plink
folder=4.3_univariate_gwas

#script=example_4.3.1_linear.sh
#echo "> $script"
#time $EXAMPLES/${folder}/${script} > /dev/null 2>&1
#rm $out_dir/*

#script=example_4.3.1_model.sh
#echo "> $script"
#time $EXAMPLES/${folder}/${script} > /dev/null 2>&1
#rm $out_dir/*

#script=example_4.3.1_assoc.sh
#echo "> $script"
#time $EXAMPLES/${folder}/${script} > /dev/null 2>&1
#rm $out_dir/*

#script=example_4.3.1_logistic.sh
#echo "> $script"
#time $EXAMPLES/${folder}/${script} > /dev/null 2>&1
#rm $out_dir/*

out_dir=
#script=example_4.3.2_lmm.sh
#echo "> $script"
#time $EXAMPLES/${folder}/${script} 
## This scripts writes output to a different directory
#rm /home/gwasuser/output/univariate_gwas/fastlm/*

# ----------------------------------------------------------------
# 4.4 
# ----------------------------------------------------------------
out_dir=/home/gwasuser/output/population_structure_correction
folder=4.4_population_structure_correction

#script=example_4.4.1.1_compute_pcs_plink.sh
#echo "> $script"
#time $EXAMPLES/${folder}/${script} > /dev/null 2>&1

#script=example_4.4.1.2_compute_pcs_eigensoft.sh
#echo "> $script"
#time $EXAMPLES/${folder}/${script} > /dev/null 2>&1

#script=example_4.4.1_correction_with_pcs.sh
#echo "> $script"
#time $EXAMPLES/${folder}/${script} > /dev/null 2>&1

#script=example_4.4.1_generate_plots.sh
#echo "> $script"
#time $EXAMPLES/${folder}/${script} > /dev/null 2>&1
#rm $out_dir/*

# ----------------------------------------------------------------
# 4.5 
# ----------------------------------------------------------------
out_dir=/home/gwasuser/output/gene_based
folder=4.5_gene_based_testing

#script=example_4.5.1.1_set_flag.sh
#echo "> $script"
#time $EXAMPLES/${folder}/${script} > /dev/null 2>&1

#script=example_4.5.1.1_make_set_flag.sh
#echo "> $script"
#time $EXAMPLES/${folder}/${script} > /dev/null 2>&1
#rm $out_dir/vegas/*

# Note: cannot run VEGAS

#script=example_4.5.2.1_fastlmm_set.sh
#echo "> $script"
#time $EXAMPLES/${folder}/${script} > /dev/null 2>&1
#rm $out_dir/plink/*

# ----------------------------------------------------------------
# 4.6 
# ----------------------------------------------------------------
out_dir=/home/gwasuser/output/epistasis
folder=4.6_epistasis

#script=example_4.6.1.1_fast_epistasis.sh
#echo "> $script"
#time $EXAMPLES/${folder}/${script} > /dev/null 2>&1

#script=example_4.6.1.2_epistasis.sh
#echo "> $script"
#time $EXAMPLES/${folder}/${script} > /dev/null 2>&1
#rm $out_dir/*

# ----------------------------------------------------------------
# 4.7 
# ----------------------------------------------------------------
out_dir=/home/gwasuser/output/network_based
folder=4.7_network_based

#script=example_4.7.1_dmgwas.sh
#echo "> $script"
#time $EXAMPLES/${folder}/${script} > /dev/null 2>&1
