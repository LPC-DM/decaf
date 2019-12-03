#!/usr/bin/env tcsh

source env_lcg.csh

pip install --user coffea

# progressbar, sliders, etc.
jupyter nbextension enable --py widgetsnbextension
