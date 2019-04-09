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

where '?' can be [1,2,3]. Fork this repo on github and git clone it:

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

The list of the inputs can be obtained in JSON format per year, packing the files beloning to the same dataset in batches of the desired size. The size of the batch will affect the performances of the following steps. The smaller the batch, the larger number of batches per dataset, the more the jobs you can run in parallel with condor. This needs to be balanced with the number of cores that can be made available through condor. A sweet spot has been found by setting the size of the batch to 75 file. 

To create the JSONS:

```
cd harvester
python pack.py --year 2017 --pack 75
```

The ```--year``` option will allow for the generation of the list for a specific year, 2017 in the example. The ```--pack``` option will set the size of the batch. The JSON files will be stored in the ```beans``` folder and they will include the cross section of the dataset each batch belongs to. A cross section of -1 is assigned to data batches.

## Generate the histograms

The generation of histograms can be launched from the ```grinder``` folder. Currently a local generation through python futures is supported together with a condor implementation.

### Running with Phyton futures

Phyton futures allows for running on multiple processor on a single node. The number of processor is currently hardcoded and set to 8. To run the local histogram generation:

```
cd ../grinder
python run.py --year 2017 --lumi 41.53 --dataset TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8____0
```

In this example, histograms for the 0 batch of the 2017 TTJets_TuneCP5_13TeV-amcatnloFXFX-pythia8 NanoAOD are being generated. The ```--year``` and the ```--lumi``` options are compulsory, ```--dataset``` is optional. Launching the script without the ```--dataset``` option will make the script run over all the batches for all the datasets. If, for example, ```--dataset TTJets``` is used, the module will run over all batches and all the datasets that match the ```TTJets``` string.

### Running with Condor

Condor will allow to parallelize jobs by running across multiple cores:

```
python submit.py --year 2017 --lumi 41.53 -t
```

This way jobs to generate the full set of 2017 histograms will be submitted to condor. the ```-t``` will allow for tarring the working environment and the necessary dependences to run on condor nodes. The module has a ```--dataset``` option that works like described before for ```run.py```. Will allow you to run on a single batch, dataset, or batches/datasets that match the input string.

