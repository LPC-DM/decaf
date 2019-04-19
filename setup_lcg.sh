#!/usr/bin/env bash

source /cvmfs/sft.cern.ch/lcg/views/LCG_94python3/x86_64-slc6-gcc62-opt/setup.sh

#pip install --user awkward
pip install --user fnal-column-analysis-tools
pip install --user awkward==0.9.0
#pip install -U --user uproot
pip install --user uproot==3.3.2
pip install --user --upgrade uproot-methods==0.3.3
# 1.14 is kindof old but pinned by other packages it seems
# pip install --user --upgrade numpy

# get dependencies for it
#pip install --user fnal-column-analysis-tools
#pip install --index-url https://test.pypi.org/simple/ --no-deps saiyan

# get latest and greatest
#git clone git@github.com:CoffeaTeam/fnal-column-analysis-tools.git
#git clone https://github.com/CoffeaTeam/fnal-column-analysis-tools.git
#git clone https://github.com/mcremone/saiyan.git

# progressbar, sliders, etc.
jupyter nbextension enable --py widgetsnbextension

# issue with python3 bindings, see https://sft.its.cern.ch/jira/browse/SPI-1198
wget https://github.com/xrootd/xrootd/archive/v4.8.4.tar.gz
tar zxf v4.8.4.tar.gz && rm v4.8.4.tar.gz
cp xrd_setup.py xrootd-4.8.4/bindings/python/
pushd xrootd-4.8.4/bindings/python/
python xrd_setup.py install --user
popd
rm -rf xrootd-4.8.4

source env_lcg.sh
