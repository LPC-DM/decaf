# Decaf

Following instructions are to run the generation of the histograms directly from NanoAOD. 

## Initial Setup

First, log into an LPC node:

```
ssh -L 9094:localhost:9094 USERNAME@cmslpc-sl7.fnal.gov
```

The command will also start forwarding the port 9094 (or whatever number you choose)to be able to use applications like jupyter once on the cluster. Then move into your area on uscms_data:

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

by running this script you will also install your grid certificate.


## Listing the input NanoAOD

The list of the inputs can be obtained in JSON format per year, packing the files beloning to the same dataset in batches of the desired size. The size of the batch will affect the performances of the following steps. The smaller the batch, the larger number of batches per dataset, the more the jobs you can run in parallel with condor. This needs to be balanced with the number of cores that can be made available through condor. A sweet spot has been found by setting the size of the batch to 75 file. 
Similarly, packing files in batches of different size will affect the performances once running with Spark. To run with Spark, use a size of 300 files per batch.

To create the JSONS:

```
cd analysis
python pack.py -y 2018 -p 300 #(or 75)
```

The ```--year``` option will allow for the generation of the list for a specific year, 2017 in the example. The ```--pack``` option will set the size of the batch. The ```--keep``` option will save .txt files that list the number of rootfiles per dataset, that will be stored in the ```metadata/YEAR``` folder. The JSON files will be stored in the ```metadata``` folder and they will include the cross section of the dataset each batch belongs to. A cross section of -1 is assigned to data batches.

## Generate the histograms

The generation of histograms can be launched from the ```analysis``` folder. Currently a local generation through python futures is supported together with a condor and a Spark implementation.

### Compile weights

Before running, weights from secondary inputs like corrections, ids ecc. need to be compiled and stored in .coffea files. To accomplish this task, simply run the following command:

```
sh compile_weights.sh 
```
the script will run the following python modules:

```
python metfilters.py
python corrections.py
python triggers.py
python ids.py
```
that of course can also run separately if only one separate set of weights needs to be compiled. Separate .coffea files will be generated, corresponding to the four python modules listed above.

### Generate Coffea Processor

The next step is generating the Coffea processor file corresponding to the analysis you want to run. For example, to generate the Dark Higgs processor, the following command can be used:

```
python generate_darkhiggs.py --year 2018
```

The module will generate a .coffea files that contains the processor instance to be used in the following steps. The ```--year``` option will allow for the generation of the processor corresponding to a specific year.

### Running with Phyton futures

Phyton futures allows for running on multiple processor on a single node. To run the local histogram generation:

```
python run.py --year 2018 --dataset MET____0_ --processor darkhiggs2018
```

In this example, histograms for the 0 batch of the 2018 MET NanoAOD are being generated. The ```--year``` option is compulsory and the ```--dataset``` is optional. The ```--lumi``` is optional. If not provided, it will default to the hard-coded lumi values in ```run.py```. Launching the script without the ```--dataset``` option will make the script run over all the batches for all the datasets. If, for example, ```--dataset TTJets``` is used, the module will run over all batches and all the datasets that match the ```TTJets``` string. The ```--processor``` option is compulsory and it will set the processor instance to be used. It corresponds to the name of a specific processor .coffea file that can be generated as described previously.

### Running with Condor

Condor will allow to parallelize jobs by running across multiple cores:

```
python submit_condor.py --year 2018 --processor darkhiggs2018 -t
```

This way jobs to generate the full set of 2018 histograms will be submitted to condor. the ```-t``` will allow for tarring the working environment and the necessary dependences to run on condor nodes. The module has a ```--dataset``` option that works like described before for ```run.py```. Will allow you to run on a single batch, dataset, or batches/datasets that match the input string.

## Running with Spark

To run with Spark, first copy your certificate to your home are on the fermicloud118 instance.

```
scp /uscms/home/USERNAME/x509up_u45169 matteoc@fermicloud118.fnal.gov:/home/USERNAME
```

Then log into fermicloud118, forwarding the port:

```
ssh -L 9094:localhost:9094 fermicloud118
```

First, check the name of the active k8s pod, where you will run your Spark applications. 

```
kubectl get pods | grep vm
```

and then, copy your certificate over the pod:

```
kubectl cp /home/USERNAME/x509up_u45169 cmsspark-vm-66-g6vbx:/home/USERNAME
```

where cmsspark-vm-66-g6vbx is the name of the pod retrieved at the previous step. Remember to port-forward to the pod too:

```
kubectl port-forward cmsspark-vm-66-g6vbx 9094:9094 &
```

Following are one-time instructions to install Laurelin, the Java library that allows for remote reading of rootfiles in Spark:

```
git clone https://github.com/lgray/laurelin.git -b useful_features
cd laurelin
rm -rf ~/.ivy2/cache/edu.vanderbilt.accre && rm -rf ~/.ivy2/jars/edu.vanderbilt.accre_laurelin*.jar && mvn install -Dmaven.test.skip=true #this nukes your jar cache and builds laurelin
kubectl cp ${HOME}/.m2/repository/edu/vanderbilt/accre/laurelin/0.4.2-SNAPSHOT/laurelin-0.4.2-SNAPSHOT.pom cmsspark-vm-66-g6vbx:${HOME}/.m2/repository/edu/vanderbilt/accre/laurelin/0.4.2-SNAPSHOT/laurelin-0.4.2-SNAPSHOT.pom
kubectl cp ${HOME}/.m2/repository/edu/vanderbilt/accre/laurelin/0.4.2-SNAPSHOT/laurelin-0.4.2-SNAPSHOT.jar cmsspark-vm-66-g6vbx:${HOME}/.m2/repository/edu/vanderbilt/accre/laurelin/0.4.2-SNAPSHOT/laurelin-0.4.2-SNAPSHOT.jar
```

To access the pod, just do:

```
/opt/ssh-into-pod.sh
```

and press 1 when asked to do so. Then:

```
su USERNAME
cd
python3 -m venv py36
source py36/bin/activate
```

If it is the first time logging into the pod, install Jupyter:

```
pip install jupyter
```
Clone decaf here:

```
git clone https://github.com/mcremone/decaf.git
```

To start Jupyter, do:

```
cd decaf
sh start_jupyter.sh
```

The script will print a link inside start_jupyter.log. Copy-paste it on your browser and you can navigate inside decaf with Jupyter. You can access the ```analysis``` folder and start the ```run_spark.ipynb``` jupyter notebook to run your analysis with Spark. At the end of your work, remember to stop Jupyter:

```
sh stop_jupyter.sh
```
