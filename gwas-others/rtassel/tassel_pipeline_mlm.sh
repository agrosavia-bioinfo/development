#!/bin/bash
# Naive
#$TASSEL_HOME/run_pipeline.pl -fork1 -plink -ped tgeno.ped -map tgeno.map \
#	-fork2 -r tpheno.tbl \
#	-combine3 -input1 -input2 -intersect -glm -export out_glm

GENOPED=tgeno.ped
GENOMAP=tgeno.map
PHENOTBL=tpheno.tbl
OUTFILE=out_mlmg

# MLM with Kinship
#$TASSEL_HOME/run_pipeline.pl \
#	-fork2 -plink -ped $GENOPED -map $GENOMAP \
#	-fork3 -r $PHENOTBL \
#    -combine4 -input2 -input3 -intersect \
#	-fork1 -plink -ped $GENOPED -map $GENOMAP -KinshipPlugin -method Dominance_Centered_IBS -endPlugin \
#	-combine5 -input4 -input1 -mlm -mlmVarCompEst P3D \
#	-mlmOutputFile $OUTFILE

# MLM with Kinship+PCs
$TASSEL_HOME/run_pipeline.pl \
	-fork1 -plink -ped $GENOPED -map $GENOMAP \
	-fork2 -r $PHENOTBL \
    -combine3 -input1 -input2 -intersect \
	-fork4 -plink -ped $GENOPED -map $GENOMAP -PrincipalComponentsPlugin -covariance -endPlugin \
	-fork5 -plink -ped $GENOPED -map $GENOMAP -KinshipPlugin -method Dominance_Centered_IBS -endPlugin \
	-combine6 -input3 -input4 -input5 -mlm -mlmVarCompEst P3D \
	-mlmOutputFile out_mlm-K+PCs
	#-combine5 -input4 -input1 -mlm -mlmVarCompEst P3D \
	#-mlmOutputFile $OUTFILE

