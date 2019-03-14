#!/usr/bin/env bash
cd pods/${1}/${3}
condor_submit run environment="year=${1}; lumi=${2}; selection=${3}; sample=${4}"
#rm run* *.tgz
cd ../../../