import cloudpickle
import pickle
import gzip
import os
from collections import defaultdict, OrderedDict
from coffea import hist, processor 
from coffea.util import load, save

def split(arr, size):
     arrs = []
     while len(arr) > size:
         pice = arr[:size]
         arrs.append(pice)
         arr   = arr[size:]
     arrs.append(arr)
     return arrs

def merge(folder,_dataset):

     filelist={}
     pd = []
     for filename in os.listdir(folder):
          if '.pkl.gz' in filename:
               if filename.split("____")[0] not in pd: pd.append(filename.split("____")[0])
    
     for pdi in pd:
          files = []
          for filename in os.listdir(folder):
               if pdi not in filename: continue
               if '.pkl.gz' not in filename: continue
               files.append(filename)
          #print(pdi,'length:',len(files))
          split_files=split(files, 100)
          #print(pdi,'number of lists:',len(split_files))
          for i in range(0,len(split_files)) :
               filelist[pdi+'___'+str(i)+'_']=split_files[i]

     coffealist=[]
     for pdi in filelist.keys():
          if _dataset not in 'None' and _dataset not in pdi: continue
          print(pdi) 
          #print(filelist[pdi])
          hists={}
          for filename in filelist[pdi]:
               fin = gzip.open(folder+'/'+filename)        
               print('Opening:',folder+'/'+filename)
               hin = cloudpickle.load(fin)
               #print('before',hin['recoil'].integrate('dataset',filename.split(".")[0]).integrate('region','isoneE').integrate('jet_selection','baggy').values())
               for k in hin.keys():
                    if k not in hists: hists[k]=hin[k]
                    else: hists[k]+=hin[k]
               #print('middle',hists['recoil'].integrate('dataset',filename.split(".")[0]).integrate('region','isoneE').integrate('jet_selection','baggy').values())
               fin.close()
               del hin
          dataset = hist.Cat("dataset", "dataset", sorting='placement')
          dataset_cats = ("dataset",)
          dataset_map = OrderedDict()
          dataset_map[pdi] = (pdi.split("___")[0]+"*",)
          for key in hists.keys():
               hists[key] = hists[key].group(dataset_cats, dataset, dataset_map)
          #print('after',hists['recoil'].integrate('dataset',pdi).integrate('region','isoneE').integrate('jet_selection','baggy').values())
          save(hists,folder+'/'+pdi+'.coffea')
          del hists
          #coffealist.append(folder+'/'+pdi+'.coffea')

     for coffeafile in os.listdir(folder):
          if '.coffea' not in coffeafile: continue
          coffealist.append(folder+'/'+coffeafile)
     print('coffealist',coffealist)

     htot={}
     for coffeafile in coffealist:
          print('Opening',coffeafile)
          hists=load(coffeafile)
          #print(hists)
          #print(coffeafile.split("/")[1].split(".")[0])
          #print('before',hists['recoil'].integrate('dataset',coffeafile.split("/")[1].split(".")[0]).integrate('region','isoneE').integrate('jet_selection','baggy').values()) 
          for k in hists:
               if k not in htot: htot[k]=hists[k]
               else: htot[k]+=hists[k]
          #print('after',htot['recoil'].integrate('dataset',coffeafile.split("/")[1].split(".")[0]).integrate('region','isoneE').integrate('jet_selection','baggy').values())
          del hists
     if _dataset in 'None': _dataset=''
     save(htot,'condor_hists_'+folder+'.coffea')

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-f', '--folder', help='folder', dest='folder')
    parser.add_option('-d', '--dataset', help='dataset', dest='dataset')
    (options, args) = parser.parse_args()

    if options.dataset: merge(options.folder,options.dataset)
    else: merge(options.folder,'None')
