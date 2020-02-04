#!/bin/bash
#/home/lg/agrosavia/opt/TASSEL5/run_pipeline.pl -fork1 -plink -ped tgeno.ped -map tgeno.map -fork2 -r tpheno.tbl -combine3 -input1 -input2 -intersect -glm -export out_glm
GENOPED=$1
GENOMAP=$2
PHENOTBL=$3
OUTFILE=$4
$TASSEL_HOME/run_pipeline.pl -fork1 -plink -ped $GENOPED -map $GENOMAP -fork2 -r $PHENOTBL -combine3 -input1 -input2 -intersect -glm -export $OUTFILE
