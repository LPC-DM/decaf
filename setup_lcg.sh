#!/usr/bin/env bash

source env_lcg.sh

pip install --user coffea==0.6.37
pip install --user https://github.com/nsmith-/rhalphalib/archive/master.zip
pip install --user xxhash
# progressbar, sliders, etc.
jupyter nbextension enable --py widgetsnbextension

