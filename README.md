# Decaf

Following instructions are to run the generation of the histograms directly from NanoAOD. 

## Initial Setup

First, log into an LPC node:

```
ssh -XY USERNAME@cmslpc-sl6.fnal.gov
```

then install your grid certificate:

```
voms-proxy-init -voms cms --valid 140:00
```

move into your area on uscms_data:

```
cd /uscms_data/d?/USERNAME/
```

where '?' can be [1,2,3]. Fork this repo on githib and git clone it:

```
git clone https://github.com/USERNAME/decaf.git
```

Then, setup the proper dependences:

```
source setup_lcg.sh
```

you will install the necessary packages as user. This is a one-time setup. When you log in next just do:

```
source env_lcg.sh
```


## Listing the input NanoAOD

