#!/usr/bin/env bash

source /cvmfs/sft.cern.ch/lcg/views/LCG_94python3/x86_64-slc6-gcc62-opt/setup.sh

#Need at least release 0.4.3 of fnal tools
pip install --user coffea

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
