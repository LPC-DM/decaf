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

def merge(folder,variable=None, exclude=None):

     lists = {}
     for filename in os.listdir(folder):
          if '.reduced' not in filename: continue
          if filename.split('--')[0] not in lists: lists[filename.split('--')[0]] = []
          lists[filename.split('--')[0]].append(folder+'/'+filename)

     for var in lists.keys():
          tmp={}
          if variable is not None:
               if not any(v==var for v in variable.split(',')): continue
          if exclude is not None:
               if any(v==var for v in exclude.split(',')): continue
          print(lists[var])
          for filename in lists[var]:
               print('Opening:',filename)
               hin = load(filename)
               if var not in tmp: tmp[var]=[hin[var]]
               else: tmp[var].append(hin[var])
               del hin
          print(tmp)
          for k in tmp:
               tmp_arr=futuresum(np.array(tmp[k]))
               hists = {}
               hists[k]=tmp_arr[0]
               print(hists)
               save(hists, folder+'/'+k+'.merged')


def postprocess(folder):
     
     variables = []
     for filename in os.listdir(folder):
          if '.merged' not in filename: continue
          #if '--' not in filename: continue
          if filename.split('.')[0] not in variables: variables.append(filename.split('.')[0])
     variables = list(filter(None, variables))

     hists = {}
     for variable in variables:
          filename = folder+'/'+variable+'.merged'
          print('Opening:',filename)
          hin = load(filename)
          hists.update(hin)
     print(hists)
     save(hists,folder+'.merged')
     
     

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-f', '--folder', help='folder', dest='folder')
    parser.add_option('-v', '--variable', help='variable', dest='variable', default=None)
    parser.add_option('-e', '--exclude', help='exclude', dest='exclude', default=None)
    parser.add_option('-p', '--postprocess', action='store_true', dest='postprocess')
    (options, args) = parser.parse_args()

    patch_mp_connection_bpo_17560()    
    if options.postprocess:
         postprocess(options.folder)
    else:
         merge(options.folder,options.variable,options.exclude)
