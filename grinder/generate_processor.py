#!/usr/bin/env python
import lz4.frame as lz4f
import cloudpickle
import pprint
import numpy as np
from coffea import hist, processor
from optparse import OptionParser
from analysis.darkhiggs import AnalysisProcessor

with lz4f.open('AnalysisProcessor.cpkl.lz4', mode='wb', compression_level=5 ) as fout:
        cloudpickle.dump(AnalysisProcessor, fout)
