# http://lcginfo.cern.ch/release/95apython3/
# Try to guess SL6 vs. CC7
if uname -r | grep -q el6; then
  source /cvmfs/sft.cern.ch/lcg/views/LCG_96python3/x86_64-slc6-gcc8-opt/setup.sh
else
  source /cvmfs/sft.cern.ch/lcg/views/LCG_96python3/x86_64-centos7-gcc8-opt/setup.sh
fi

export PYTHONPATH=~/.local/lib/python3.6/site-packages:$PYTHONPATH
export PYTHONPATH=$PWD/analysis:$PYTHONPATH
export PYTHONWARNINGS="ignore"
#export PATH=~/.local/bin:$PATH
cat $HOME/private/$USER.txt | voms-proxy-init -voms cms --valid 140:00 -pwstdin
