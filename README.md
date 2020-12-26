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
cd analysis/metadata
python pack.py -y 2018 -p 300 #(or 75)
```

The ```--year``` option will allow for the generation of the list for a specific year, 2017 in the example. The ```--pack``` option will set the size of the batch. The ```--keep``` option will save .txt files that list the number of rootfiles per dataset, that will be stored in the ```metadata/YEAR``` folder. The JSON files will be stored in the ```metadata``` folder and they will include the cross section of the dataset each batch belongs to. A cross section of -1 is assigned to data batches.

## Generate the histograms

The generation of histograms can be launched from the ```analysis``` folder. Currently a local generation through python futures is supported together with a condor and a Spark implementation.

### Compile weights

Before running, weights from secondary inputs like corrections, ids ecc. need to be compiled and stored in .coffea files. To accomplish this task, simply run the following command:

```
sh generate_secondary_inputs.sh
```
the script will run the following python modules:

```
python secondary_inputs/metfilters.py
python secondary_inputs/corrections.py
python secondary_inputs/triggers.py
python secondary_inputs/ids.py
```
these modules can also run separately by passing to the bash script the argument corresponding to the name of the module you want to run. For example, to run ```secondary_inputs/corrections.py```, just do:

```
sh generate_secondary_inputs.sh corrections
```

Separate .coffea files will be generated, corresponding to the four python modules listed above, and stored in the ```secondary_inputs``` folder.

### Generate Coffea Processor

The next step is generating the Coffea processor file corresponding to the analysis you want to run. For example, to generate the Dark Higgs processor, the following command can be used:

```
sh generate_processor.sh darkhiggs 2018
```

The script will run the following command, using the first argument to select the corresponding python module in the ```processors``` folder and the second argument to set the ```--year``` option:

```
python processors/darkhiggs.py --year 2018
```

The module will generate a .coffea file, stored in the ```processors``` folder, that contains the processor instance to be used in the following steps. The ```--year``` option will allow for the generation of the processor corresponding to a specific year.

### Running with Python futures

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
python3 -m venv py36 --system-site-packages
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

Before starting Jupyter, set the environmental proxy variable to your certificate:

```
export X509_USER_PROXY=/home/matteoc/x509up_u45169
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

## From Coffea Histograms to Fits

### Generating the model

Taking as example the dark Higgs analysis, run the following command to generate the background model:

```
python models/darkhiggs.py -y 2018
```
The ```models/darkhiggs.py``` module extracts coffea histgrams from the ```hists/darkhiggs201?.scaled``` files and utilize them to generate the different templates that will later be rendered into datacards/workspaces. It also defines transfer factors for the data-driven background models. It produces ```.model``` files that are saved into the ```data``` folder. More specifically, the ```models/darkhiggs.py``` module produces one ```.model``` file for each pass/fail and recoil bin. 

### Rendering the model

The models saved in the ```.model``` files at the previous step can be rendered into datacards by running the following command:

```
python render.py -m model_name
```

the ```-m``` or ```--model``` options provide in input the name of the model to render, that corresponds to the name of the ```.model``` file where the model is stored. The ```render.py``` module launches python futures jobs to process in parallel the different control/signal regions that belongs to a single model. Different models can also be rendered in parallel, by using condor. In order to run rendering condor job, the following command should be ran:

```
python render_condor.py -m model_name
```

this time, all the models stored into ```.model``` files whose name contains the sting passed via the  ```-m``` options are going to be rendered.
In order for the rendering step to be completed successfully and for the resulting datacards and workspaces to be used in later steps, both the python version and the ROOT version should correspond to the ones that are going to be set when running ```combine```. It is therefore suggested to clone the ```decaf``` repository inside the ```src``` folder of ```CMSSSW``` release that will be used to perform the combine fit and run the rendering from there, without running the ```env_lcg.sh```. You should, instead, do the regular ```cmsenv``` command.
At the time of writing, the recommended version is ```CMSSW_10_2_13``` (https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/#cc7-release-cmssw_10_2_x-recommended-version). This version uses python 2.7 and ROOT 6.12.
In addition, you should install the following packages (to be done *after* ```cmsenv```):

```
pip install --user flake8
pip install --user cloudpickle
pip install --user https://github.com/nsmith-/rhalphalib/archive/master.zip
```

Finally, since we are using Python2.7, you should edit the file
```
$HOME/.local/lib/python2.7/site-packages/rhalphalib/model.py
```
and change line 94 from
```
obsmap.insert(ROOT.std.pair('const string, RooDataHist*')(channel.name, obs))
```
to
```
obsmap.insert(ROOT.std.pair('const string, RooDataHist*')(channel.name.encode("ascii"), obs))
```

The ```render.py``` module saves datacards and workspaces into folders that show the following naming:

```
datacards/model_name
```

The ```render_condor.py``` module returns ```.tgz``` tarballs that contain the different datacards/workspaces, and are stored into the ```datacards``` folder. To untar them, simply do:

```
python macros/untar_cards.py -a darkhiggs
```

Where the ```-a``` or ```--analysis``` options correspond to the analysis name. The ```untar_cards.py``` script will untar all the tarballs that contain the string that is passed through the ```-a``` option.
To merge the datacards, run the following script:

```
python macros/combine_datacards.py -a darkhiggs
```
Where the ```-a``` or ```--analysis``` options correspond to the analysis name. The ```combine_datacards.py``` script will combine all the datacards whose name contains the string that is passed through the ```-a``` option. The script will create a folder inside ```datacards``` whose name corresponds to the string that is passed through the ```-a``` option, will move all the workspaces that correspond to the datacards it combined inside it, and will save in it the combined datacard, whose name will be set to the string that is passed through the ```-a``` option.

### Using Combine

Move inside the newly generated folder:

```
cd datacards/darkhiggs
```

From here, to convert the datacard into the workspace you will use to run the fit, do:

```
text2workspace.py darkhiggs.txt
```

The result of this step will be a ```darkhiggs.root``` file to be used for the combine fit.
where ```darkhiggs.txt``` is the name of the combined datacard generated at the previous step.

#### Running special datacards

Before running the combie tools, you have to put ```ulimit -s unlimited``` command.

When you render the model, you can use the alternative datacards:
```
analysis/models/darkhiggs_fourregions.py
analysis/models/darkhiggs_fiveregions.py
analysis/models/darkhiggs_eightregions.py
```

that try the fit with "Signal, WMuon, TopMuon, DoubleMuon" (four);
with all those plus "Gamma" (five);
with all those plus "WEle, TopEle, DoubleEle" (eight).
In all of those the minor MCs are not added at all to the fit.
The usual `darkhiggs.py` has all eight regions + minor MCs.

Alternative fits:

Minuit:
```
combine -M FitDiagnostics -d darkhiggs.root \
--expectSignal 1 --forceRecreateNLL --cminDefaultMinimizerType Minuit \
--ignoreCovWarning --saveShapes --saveWithUncertainties --saveWorkspace --plots \
&> out_FitDiagnostics_Minuit.log
```

Minuit (fast):   
Drop the options: `--saveShapes`, `--saveWorkspace`, `--plots`   
```
combine -M FitDiagnostics -d darkhiggs.root \
--expectSignal 1 --forceRecreateNLL --cminDefaultMinimizerType Minuit \
--ignoreCovWarning --saveWithUncertainties &> out_FitDiagnostics_Minuit.log
```

Robust fit:
```
combine -M FitDiagnostics -d darkhiggs.root \
--expectSignal 1 --forceRecreateNLL  --robustFit=1 \
--ignoreCovWarning --saveShapes --saveWithUncertainties --saveWorkspace --plots \
&> out_FitDiagnostics_robustFit.log
```

## Macros

### Usage examples

```
python macros/dump_templates.py -w datacards/darkhiggs/darkhiggs.root:w --observable fjmass -o plots/darkhiggs2018/model
```
```
python macros/hessian.py -w datacards/darkhiggs/higgsCombineTest.FitDiagnostics.mH120.root:w -f datacards/darkhiggs/fitDiagnostics.root:fit_b
```
#### Pulls plotting 
Change the path to store outputs in `plotConfig.py`    
```
python macros/diffNuisances.py -g pulls.root ../fitDiagnostics.root
```

#### Postfit plotting 
```
python macros/plotStackedPostFit_darkhiggs.py 
```
