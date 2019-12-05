#!/usr/bin/env bash

source env_lcg.sh

pip install --user coffea
pip install --user https://github.com/nsmith-/rhalphalib/archive/master.zip
# progressbar, sliders, etc.
jupyter nbextension enable --py widgetsnbextension

