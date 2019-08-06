# http://lcginfo.cern.ch/release/94python3/
unsetenv PYTHON_LIB_SITE_PACKAGES
set lsalias=`which ls | sed 's|.*to \(.*\)|\1|'`
unalias ls
source /cvmfs/sft.cern.ch/lcg/views/LCG_94python3/x86_64-slc6-gcc62-opt/setup.csh
alias ls $lsalias

#export PATH=${HOME}/.local/bin:$PATH
#export PYTHONPATH=${HOME}/.local/lib/python3.6/site-packages:$PYTHONPATH
