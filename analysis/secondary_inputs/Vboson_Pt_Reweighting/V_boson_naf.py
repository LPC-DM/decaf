import sys
import ROOT
import stat
import os
import Utilities.General.cmssw_das_client as das_client

file_prefix = "root://xrootd-cms.infn.it//"

veto_list = [
    # 2018
    "/ZJetsToNuNu_HT-400To600_TuneCP5_13TeV-madgraph/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v1/MINIAODSIM",
    "/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v2/MINIAODSIM",
    "/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v3/MINIAODSIM",
    "/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15-v4/MINIAODSIM",
    "/DYJetsToLL_M-50_HT-400to600_TuneCP5_PSweights_13TeV-madgraphMLM-pythia8/RunIIAutumn18MiniAOD-102X_upgrade2018_realistic_v15_ext2-v3/MINIAODSIM",
    # 2017
    "/ZJetsToNuNu_HT-1200To2500_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    "/ZJetsToNuNu_HT-400To600_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    "/ZJetsToNuNu_HT-600To800_13TeV-madgraph/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    "/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    "/DYJetsToLL_M-50_HT-100to200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",
    "/DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    "/DYJetsToLL_M-50_HT-200to400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",
    "/DYJetsToLL_M-50_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    "/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    "/DYJetsToLL_M-50_HT-400to600_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",
    "/DYJetsToLL_M-50_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    "/DYJetsToLL_M-50_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    "/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_1core_94X_mc2017_realistic_v14_ext1-v1/MINIAODSIM",
    "/GJets_HT-40To100_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v2/MINIAODSIM",
    "/GJets_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
    "/GJets_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8/RunIIFall17MiniAODv2-PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/MINIAODSIM",
]


def get_files(dataset_name):
    print (dataset_name)
    data = das_client.get_data("dataset=" + dataset_name + " instance=prod/global").get(
        "data", None
    )
    datasets = [
        data[i].get("dataset", None)[0].get("name", None) for i in range(len(data))
    ]
    print (datasets)
    files = []
    for dataset in datasets:
        if dataset in veto_list:
            continue
        print (dataset)
        data = das_client.get_data("file dataset=" + dataset + " instance=prod/global")
        for d in data.get("data", None):
            # print(d)
            for f in d.get("file", None):
                # print(f)
                # if not 'nevents' in f:
                # continue
                files.append(
                    [file_prefix + f.get("name", None), f.get("nevents", None)]
                )
    return files


def split_files_into_jobs(files, events_per_job):
    events = 0
    file_splitting = []
    files_in_job = []
    for i, file in enumerate(files):
        # print("file ",i)
        files_in_job.append(file[0])
        if file[1] == None:
            file_ = ROOT.TFile.Open(file)
            tree = None
            try:
                tree = file_.Get("Events")
            except ReferenceError:
                file_.Close()
                continue
            events += tree.GetEntries()
            file_.Close()
        else:
            events += file[1]
        if events >= events_per_job or i == (len(files) - 1):
            file_splitting.append(files_in_job)
            events = 0
            files_in_job = []
    return file_splitting


def print_shell_script(boson, postfix, files, era):
    script = ""
    script += "#!/bin/bash\n"
    script += "export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch\n"
    script += "source $VO_CMS_SW_DIR/cmsset_default.sh\n"
    script += "cd /nfs/dust/cms/user/mwassmer/MonoTop/CMSSW_10_2_18/src\n"
    script += "eval `scram runtime -sh`\n"
    if not os.path.isdir("root_files"):
        os.mkdir("root_files")
    script += "cd /nfs/dust/cms/user/mwassmer/MonoTop/useful_scripts/Vboson_Pt_Reweighting/root_files\n"
    script += (
        "python /nfs/dust/cms/user/mwassmer/MonoTop/useful_scripts/Vboson_Pt_Reweighting/V_boson_pt_reweighting.py "
        + era
        + " "
        + boson
        + " "
        + postfix
    )
    for file in files:
        script += " " + file
    script += "\n"
    script += "exitcode=$?\n"
    script += "#" + boson + "_boson_pt_" + postfix + ".sh\n"
    script += "if [ $exitcode -eq 0 ]\n"
    script += "then\n"
    script += "  exit 0\n"
    script += "else\n"
    script += "  exit 1\n"
    script += "fi\n"

    if not os.path.isdir("scripts"):
        os.mkdir("scripts")
    filename = "scripts/" + boson + "_boson_pt_" + era + "_" + postfix + ".sh"
    f = open(filename, "w")
    f.write(script)
    f.close()
    print ("created script", filename)
    st = os.stat(filename)
    os.chmod(filename, st.st_mode | stat.S_IEXEC)


era = str(sys.argv[1])
boson = str(sys.argv[2])

if not (boson == "Zvv" or boson == "Zll" or boson == "W" or boson == "G"):
    print ("first argument has to be Zll or Zvv or W or G (photon)")
    exit()

files = get_files(str(sys.argv[3]).replace('"', ""))
print ("number of files: ", len(files))
print ("number of events: ", sum(element[1] for element in files))
file_splitting = split_files_into_jobs(files, 1000000)
# print(file_splitting)
for i, files in enumerate(file_splitting):
    print_shell_script(boson, str(i), files, era)
