#!/bin/bash
#
# Create multiGWAS profile according to current directory
# 
#
echo "Creating multiGWAS profile..."

echo ""
echo "#------------------- multiGWAS HOME ---------------------" >  multiGWAS-profile.sh
cmm1="export MULTIGWAS_HOME=$PWD" 
echo $cmm1 >> multiGWAS-profile.sh; eval $cmm1
cmm2='MULTIGWAS_TOOLS=$MULTIGWAS_HOME/tools' 
echo $cmm2 >> multiGWAS-profile.sh; eval $cmm2
cmm3='MULTIGWAS_SOURCES=$MULTIGWAS_HOME/sources'
echo $cmm3 >> multiGWAS-profile.sh; eval $cmm3
cmm4='export PATH=$PATH:$MULTIGWAS_TOOLS:$MULTIGWAS_SOURCES'
echo $cmm4 >> multiGWAS-profile.sh; eval $cmm4


echo ""
echo "#------------------- Java config ---------------------" >>  multiGWAS-profile.sh
echo "export JAVA_HOME=$PWD/opt/java" >> multiGWAS-profile.sh
echo 'export PATH=$PATH:$JAVA_HOME/bin'>> multiGWAS-profile.sh
export JAVA_HOME=$PWD/opt/java
export PATH=$PATH:$JAVA_HOME/bin

echo ""
echo "#------------------- R config ---------------------" >>  multiGWAS-profile.sh
echo "R_HOME=$PWD/opt/R" >> multiGWAS-profile.sh
echo 'export PATH=$R_HOME/bin:$PATH'>> multiGWAS-profile.sh
R_HOME=$PWD/opt/R
 
echo ""
echo "#------------------- multiGWAS.R tool profile ---------------------" >> $HOME/.profile
echo "source $PWD/multiGWAS-profile.sh" >> $HOME/.profile


#------------------- Ubuntu dependencies: r-cran-devtools and libcurl4-openssl-dev ---------------------" 
sudo apt update
sudo apt-get install build-essential libcurl4-openssl-dev gfortran zlib1g-dev libbz2-dev liblzma-dev libpcre3-dev add-apt-key software-properties-common

echo ""
#------------------- R install ---------------------" 
cd $PWD/opt/installers/R-3.6.1

echo ""; echo ""; echo ""
read -p ">>>>>>>>>>>>>>> Press any key to configure/install locally R 6.1 (Ctrl-C to abort)" 
echo ""; echo ""; echo ""
./configure --prefix=$R_HOME --enable-R-shlib --enable-java=no --with-x=no --with-cairo --with-readline=no
echo ""; echo ""; echo ""
read -p ">>>>>>>>>>>>>>> Press any key to compile (make) R 6.1 (Ctrl-C to abort)" 
echo ""; echo ""; echo ""
make

echo ""; echo ""; echo ""
read -p ">>>>>>>>>>>>>>> Press any key to install (make install) R 6.1 (Ctrl-C to abort)" 
echo ""; echo ""; echo ""
make install

echo "#------------------- multiGWAS.R tool profile ---------------------" >> $HOME/.profile
echo "Updating user profile to execute multiGWAS-profile at startup..."
echo "source $MULTIGWAS_HOME/multiGWAS-profile.sh" >> $HOME/.profile

echo ""
echo "MultiGWAS is ready to use, right after installed!"
echo ""
echo ""


# OTHER elements to install
# For PNGs and VennDiagram transparency
# sudo apt-get install libcairo2-dev
# Buid r using: "configure --prefix=/home/mg/multiGWAS/opt/R --with-x=yes --with-cairo"

# For rmarkdown
# sudo apt-get install pandoc

# For markdown tables
# sudo apt-get install libxml2-dev
# sudo apt-get install libssl-dev

