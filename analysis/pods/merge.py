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

def merge(directory):

     filelist={}
     pd = []
     for filename in os.listdir(directory):
          if '.pkl.gz' in filename:
               if filename.split("____")[0] not in pd: pd.append(filename.split("____")[0])
    
     for pdi in pd:
          files = []
          for filename in os.listdir(directory):
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
          print(pdi) 
          #print(filelist[pdi])
          hists={}
          for filename in filelist[pdi]:
               fin = gzip.open(directory+'/'+filename)        
               print('Opening:',directory+'/'+filename)
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
          save(hists,directory+'/'+pdi+'.coffea')
          del hists
          coffealist.append(directory+'/'+pdi+'.coffea')
          
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
     save(htot,'condor_hists'+directory+'.coffea')

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-d', '--directory', help='directory', dest='directory')
    (options, args) = parser.parse_args()

    merge(options.directory)
