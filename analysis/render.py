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
        small_model = rl.Model('darkhiggs_'+str(ch.name))
        small_model.addChannel(model[str(ch.name)])
        model_arr.append(small_model)
    print(model_arr)
    print('Rendering')

    with concurrent.futures.ProcessPoolExecutor(max_workers=len(model_arr)) as executor:
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
    os.system('mkdir -p datacards/'+options.model)
    render(options.model)
