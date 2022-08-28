import concurrent.futures
import cloudpickle
import pickle
import gzip
import os
import numpy as np
from collections import defaultdict, OrderedDict
from coffea import hist, processor 
from coffea.util import load, save
from helpers.futures_patch import patch_mp_connection_bpo_17560

def add(chunk_tmp_arr):
     print('Job started')
     return np.sum(chunk_tmp_arr)

def futuresum(tmp_arr):
     while np.size(tmp_arr)>1:
          chunk_sum=[]
          chunk_tmp_arr = np.array_split(tmp_arr, int(np.size(tmp_arr)/2))
          if len(chunk_tmp_arr)>1:
               with concurrent.futures.ProcessPoolExecutor(max_workers=16) as executor:
                    futures = set()
                    futures.update(executor.submit(add,chunk_tmp_arr[i]) for i in range(0,len(chunk_tmp_arr)))
                    if(len(futures)==0): continue
                    try:
                         total = len(futures)
                         processed = 0
                         while len(futures) > 0:
                              finished = set(job for job in futures if job.done())
                              for job in finished:
                                   chunk_i = job.result()
                                   chunk_sum.append(chunk_i)
                              futures -= finished
                         del finished
                    except KeyboardInterrupt:
                         print("Ok quitter")
                         for job in futures: job.cancel()
                    except:
                         for job in futures: job.cancel()
                         raise
          else:
               chunk_sum.append(add(chunk_tmp_arr[0]))
          tmp_arr=np.array(chunk_sum)
          print(tmp_arr)
     return tmp_arr


def reduce(folder,_dataset=None,variable=None):

     lists = {}
     for filename in os.listdir(folder):
          if '.futures' not in filename: continue
          if filename.split("____")[0] not in lists: lists[filename.split("____")[0]] = []
          lists[filename.split("____")[0]].append(folder+'/'+filename)
          
     for pdi in lists.keys():
          if _dataset is not None and _dataset not in pdi: continue
          tmp={}
          for filename in lists[pdi]:
               print('Opening:',filename)
               hin = load(filename)
               for k in hin.keys():
                    if variable is not None and k!=variable: continue
                    print('Considering variable',k)
                    if k not in tmp: tmp[k]=[hin[k]]
                    else: tmp[k].append(hin[k])
               del hin
          for k in tmp:
               tmp_arr=futuresum(np.array(tmp[k]))
               hists = {}
               hists[k]=tmp_arr[0]
               dataset = hist.Cat("dataset", "dataset", sorting='placement')
               dataset_cats = ("dataset",)
               dataset_map = OrderedDict()
               for d in hists[k].identifiers('dataset'):
                    if d.name.split("____")[0] not in dataset_map: dataset_map[d.name.split("____")[0]] = (d.name.split("____")[0]+"*",)
               hists[k] = hists[k].group(dataset_cats, dataset, dataset_map)
               print(hists)
               save(hists, folder+'/'+k+'--'+pdi+'.reduced')

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-f', '--folder', help='folder', dest='folder')
    parser.add_option('-d', '--dataset', help='dataset', dest='dataset', default=None)
    parser.add_option('-v', '--variable', help='variable', dest='variable', default=None)
    (options, args) = parser.parse_args()

    patch_mp_connection_bpo_17560()    
    reduce(options.folder,options.dataset,options.variable)
