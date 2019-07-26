#!/usr/bin/env bash

source env_lcg.sh

pip install --user coffea
pip -q install https://github.com/Parsl/parsl/zipball/master

# progressbar, sliders, etc.
jupyter nbextension enable --py widgetsnbextension

