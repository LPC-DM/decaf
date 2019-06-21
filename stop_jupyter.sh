#!/bin/bash
jupyter notebook stop 9094
condor_rm matteoc -name lpcschedd1.fnal.gov
condor_rm matteoc -name lpcschedd2.fnal.gov
condor_rm matteoc -name lpcschedd3.fnal.gov
