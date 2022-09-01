from __future__ import print_function, division
from optparse import OptionParser
from collections import defaultdict, OrderedDict
import concurrent.futures
import sys
import os
import rhalphalib as rl
import gzip
import pickle
import cloudpickle
import ROOT
import gzip
sys.setrecursionlimit(1500)
rl.util.install_roofit_helpers()
rl.ParametericSample.PreferRooParametricHist = False

def futurerender(some_model, directory):
    print('Job started for',some_model)
    return some_model.renderCombine('datacards/'+directory)

def render(modelname):
    
    with open('data/'+modelname+'.model') as fin:
        model = pickle.load(fin)

    model_arr = []
    for ch in model:
        print('generating model for channel',ch.name)
        small_model = rl.Model(ch.name.encode("ascii"))
        small_model.addChannel(model[ch.name])
        model_arr.append(small_model)
    print(model_arr)
    print('Rendering')

    #for i in range(0,len(model_arr)):
    #    futurerender(model_arr[i], modelname)

    with concurrent.futures.ProcessPoolExecutor(max_workers=16) as executor:
        futures = set()
        futures.update(executor.submit(futurerender,model_arr[i], modelname) for i in range(0,len(model_arr)))
        try:
            total = len(futures)
            processed = 0
            while len(futures) > 0:
                finished = set(job for job in futures if job.done())
                for job in finished:
                    job.result()
                futures -= finished
            del finished
        except KeyboardInterrupt:
            print("Ok quitter")
            for job in futures: job.cancel()
        except:
            for job in futures: job.cancel()
            raise

if __name__ == '__main__':
    if not os.path.exists('datacards'):
        os.mkdir('datacards')
    parser = OptionParser()
    parser.add_option('-m', '--model', help='model', dest='model', default='')
    (options, args) = parser.parse_args()
    for filename in os.listdir('data'):
        if '.model' not in filename: continue
        if options.model:
            if not any(model in filename for model in options.model.split(',')): continue
        os.system('mkdir -p datacards/'+filename.split('.')[0])
        render(filename.split('.')[0])
