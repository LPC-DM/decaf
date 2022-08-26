#!/bin/bash

echo "Make model 2016"
python models/darkhiggs.py -y 2016 -f -m 40to120
python models/darkhiggs.py -y 2016 -f -m 120to300
echo ""

echo "Make model 2017"
python models/darkhiggs.py -y 2017 -f -m 40to120
python models/darkhiggs.py -y 2017 -f -m 120to300
echo ""

echo "Make model 2018"
python models/darkhiggs.py -y 2018 -f -m 40to120
python models/darkhiggs.py -y 2018 -f -m 120to300
