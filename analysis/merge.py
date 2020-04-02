import concurrent.futures
import cloudpickle
import pickle
import gzip
import os
from collections import defaultdict, OrderedDict
from coffea import hist, processor 
from coffea.util import load, save

'''
def patch_future(cls):
     def __iter__(self):
          if not self.done(): 
               yield self
          return self.result()
     cls.__iter__ = __iter__

from concurrent.futures import Future
patch_future(Future)
'''

def split(arr, size):
     arrs = []
     while len(arr) > size:
         pice = arr[:size]
         arrs.append(pice)
         arr   = arr[size:]
     arrs.append(arr)
     return arrs

def merging(mergingfile,filelist):
     print('Generating merged file:',mergingfile+'.merged') 
     hists={}
     for filename in filelist:
          print('Opening:',filename)
          hin = load(filename)
          for k in hin.keys():
               if k not in hists: hists[k]=hin[k]
               else: hists[k]+=hin[k]
          #fin.close()
          del hin
     dataset = hist.Cat("dataset", "dataset", sorting='placement')
     dataset_cats = ("dataset",)
     dataset_map = OrderedDict()
     for d in hists[list(hists.keys())[0]].identifiers('dataset'):
          new_dname = d.name.split("____")[0] + '____' + mergingfile.split("____")[1]
          print("Merging",d.name,"into",new_dname)
          if new_dname not in dataset_map:
               dataset_map[new_dname] = (d.name.split("____")[0]+"*",)
     for key in hists.keys():
          hists[key] = hists[key].group(dataset_cats, dataset, dataset_map)
     return mergingfile, hists

def merge(folder,_dataset=None):

     pd = []
     for filename in os.listdir(folder):
          if '.futures' in filename:
               if filename.split("____")[0] not in pd: pd.append(filename.split("____")[0])
     print('List of primary datasets:',pd)

     for pdi in pd:
          if _dataset is not None and _dataset not in pdi: continue
          files = []
          for filename in os.listdir(folder):
               if '.futures' not in filename: continue
               if pdi not in filename: continue
               files.append(folder+'/'+filename)
          split_files=split(files, int(len(files)/8) )
          with concurrent.futures.ProcessPoolExecutor(max_workers=len(split_files)) as executor:
               futures = set()
               futures.update(executor.submit(merging,pdi+'____'+str(i)+'_',split_files[i]) for i in range(0,len(split_files)) )
               if(len(futures)==0): continue
               try:
                    total = len(futures)
                    processed = 0
                    while len(futures) > 0:
                         finished = set(job for job in futures if job.done())
                         for job in finished:
                              mergingfile, hists = job.result()
                              save(hists,folder+'/'+mergingfile+'.merged')
                              del hists
                              print("Processing: done with % 4d / % 4d files" % (processed, total))
                              processed += 1
                         futures -= finished
                    del finished
               except KeyboardInterrupt:
                    print("Ok quitter")
                    for job in futures: job.cancel()
               except:
                    for job in futures: job.cancel()
                    raise

     mergedlist=[]
     for mergedfile in os.listdir(folder):
          if '.merged' not in mergedfile: continue
          mergedlist.append(folder+'/'+mergedfile)
     print('List of merged files:',mergedlist)

     htot={}
     for mergedfile in mergedlist:
          print('Opening:',mergedfile)
          hists=load(mergedfile)
          for k in hists:
               if k not in htot: htot[k]=hists[k]
               else: htot[k]+=hists[k]
          del hists
     save(htot,folder+'.merged')

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-f', '--folder', help='folder', dest='folder')
    parser.add_option('-d', '--dataset', help='dataset', dest='dataset')
    (options, args) = parser.parse_args()

    dataset=None
    if options.dataset: dataset=options.dataset
    merge(options.folder,dataset)
