#!/bin/bash
if [ -z "$@" ]
then
    python secondary_inputs/corrections.py
    python secondary_inputs/metfilters.py
    python secondary_inputs/triggers.py
    python secondary_inputs/ids.py
else
    python secondary_inputs/$1.py
fi
