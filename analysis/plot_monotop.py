from __future__ import print_function, division
from collections import defaultdict, OrderedDict
import gzip
import pickle
import json
import os
import uproot
import matplotlib.pyplot as plt

import numpy as np
from coffea import hist, processor
from coffea.hist import plot
import pdb

hists={}
pd = []
year = '2018'
dirname = 'pods/monotop' + year
for filename in os.listdir(dirname):
    #if 'MET' in filename or 'SingleElectron' in filename or 'SinglePhoton' in filename or 'EGamma' in filename: continue
    if '.pkl.gz' in filename:
        if filename.split("____")[0] not in pd: pd.append(filename.split("____")[0])
        with gzip.open(dirname+'/'+filename) as fin:
            hin = pickle.load(fin)
            for k in hin.keys():
                #if 'e1drj' in k or 'sumw' in k or 'recoil' in k:
                if k in hists: hists[k]+=hin[k]
                else: hists[k]=hin[k]

pdataset = hist.Cat("pdataset", "pdataset", sorting='placement')
pdataset_cats = ("dataset",)
pdataset_map = OrderedDict()

for pdi in pd:
    pdataset_map[pdi] = (pdi+"*",)

for key in hists.keys():
    hists[key] = hists[key].group(pdataset_cats, pdataset, pdataset_map)

scale={}
for pdi in hists['sumw'].identifiers('pdataset'):
    scale[pdi]=hists['sumw'].project('pdataset',overflow='all').values()[(str(pdi),)]

for key in hists.keys():
    if key=='sumw': continue
    for pdi in hists[key].identifiers('pdataset'):
        hists[key].scale({pdi:1/scale[pdi]},axis='pdataset')

data_hists={}

for filename in os.listdir(dirname):
    #    if 'MET' in filename or 'SingleElectron' in filename or 'SinglePhoton' in filename or 'EGamma' in filename:
    if '.pkl.gz' in filename:
        with gzip.open(dirname+'/'+filename) as fin:
            hin = pickle.load(fin)
            for k in hin.keys():
                #if hin[k].identifiers('category')[0] not in 'isoneM': continue
                if k in data_hists: data_hists[k]+=hin[k]
                else: data_hists[k]=hin[k]

process = hist.Cat("process", "Process", sorting='placement')
process_cats = ("pdataset",)
process_map = OrderedDict()

process_map["WW"] = ("WW*",)
process_map["WZ"] = ("WZ*",)
process_map["ZZ"] = ("ZZ*",)
#process_map["Hbb"] = ("*HToBB*")
process_map["QCD"] = ("QCD*",)
process_map["ST"] = ("ST*",)
process_map["TT"] = ("TTJets*",)
process_map["Wjets"] = ("WJets*",)
process_map["Zjets"] = ("DYJets*",)

#process_map["Znunu"] = ("ZJets*",)

data_r_map = OrderedDict()
data_r_map['isoneE'] = 'EGamma'
data_r_map['isoneM'] = 'MET'
data_r_map['istwoE'] = 'EGamma'
data_r_map['istwoM'] = 'MET'
#data_r_map['isoneA'] = 'SinglePhoton'
data_r_map['hadronic'] = 'MET'

data_map = OrderedDict()
data_map["MET"] = ("MET*", )
data_map["EGamma"] = ("EGamma*", )
#data_map["SinglePhoton"] = ("EGamma*", )
data_cats = ("dataset",)
print("Plotting Histograms")
for key in hists.keys():
    hists[key] = hists[key].group(process_cats, process, process_map)
    data_hists[key] = data_hists[key].group(data_cats, process, data_map)

hists['recoil'].axis('recoil').label = 'Hadronic Recoil (GeV)'
hists['fj1pt'].axis('fj1pt').label = 'AK15 Leading Jet Pt (GeV)'
hists['j1pt'].axis('j1pt').label = 'AK4 Leading Jet Pt (GeV)'
hists['fjmass'].axis('fjmass').label = 'AK15 Leading Jet Mass (GeV)'

#print(hists['recoil'].project('process','Hbb').values())

#data_map['isoneE'] = 'SingleElectron'
data_map['isoneE'] = 'EGamma'
data_map['isoneM'] = 'MET'
data_map['istwoE'] = 'EGamma'
data_map['istwoM'] = 'MET'
#data_map['isoneA'] = 'SinglePhoton'
#data_map['iszeroL'] = 'MET'
data_map['hadronic'] = 'MET'

print(hists['recoil'].identifiers('category'))
for r in hists['recoil'].identifiers('category'):
    exp = 0
    print('------------------')
    print('------------------')
    print('Category:',r)
    print('------------------')
    for p in hists['recoil'].integrate('category',r).identifiers('process'):
        for s in hists['recoil'].integrate('category',r).integrate('process',p).identifiers('control_region'):
            #pdb.set_trace()
            yld = np.sum(hists['recoil'].values(overflow='all')[(str(p),str(r),str(s))])
            exp += yld
            print('Category',r,'Process:',p,'Control Region:',s, '%.1f' % yld)
    print('------------------')
    print('Total expected:', '%.1f' % exp)
    print(data_hists['recoil'].integrate('category',r).identifiers('control_region'))
    for s in data_hists['recoil'].integrate('process',data_map[str(r)]).integrate('category',r).identifiers('control_region'):
        print('Category',r,'Dataset:',data_map[str(r)],'Control Region:',s)
        #pdb.set_trace()
        print('Total observed(',data_map[str(r)],str(r),str(s),'):', '%.1f' % np.sum(data_hists['recoil'].values(overflow='all')[(data_map[str(r)],str(r),str(s))]))
    else:
        print("Category",r,"not found in data_hists category")
    print('------------------')
    print('------------------')
    print()

from cycler import cycler

plt.rcParams.update({'font.size': 14, 'axes.titlesize': 18, 'axes.labelsize': 18, 'xtick.labelsize': 12, 'ytick.labelsize': 12})
fill_opts = {'edgecolor': (0,0,0,0.3), 'alpha': 0.8}
error_opts = {'label':'Stat. Unc.', 'hatch':'///', 'facecolor':'none', 'edgecolor':(0,0,0,.5), 'linewidth': 0}
nostack_fill_opts = {'alpha': 0.2, 'label': '_nolabel_'}
#data_err_opts = {'linestyle':'none', 'marker': '.', 'markersize': 10., 'color':'k', 'elinewidth': 1, 'emarker': '_'}
data_err_opts = {'linestyle':'none', 'marker': '.', 'markersize': 10., 'color':'k', 'elinewidth': 1}
colors = ['#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',"#ff8000","#7f00ff","#000000"]
cat_names = {'isoneM':'Single Muon', 'isoneE':'Single Electron', 'istwoM':'Double Muon', 'istwoE':'Double Electron', 'hadronic':'All Hadronic'}

if not os.path.exists('stack'): os.makedirs('stack')
for key in hists.keys():
    if key=='sumw': continue

    for cat in hists[key].identifiers('category'):
        for s in hists[key].integrate('category',cat).identifiers('control_region'):
            #pdb.set_trace()
            print("hist:",key,"cat:",cat,"CR:",s)
            if str(s)=="QCD": continue
            fig, (ax, rax) = plt.subplots(2, 1, figsize=(11,13), gridspec_kw={"height_ratios": (3, 1)}, sharex=True)
            fig.subplots_adjust(hspace=.07)
            ax.set_prop_cycle(cycler(color=colors))
            if hists[key].integrate('category',cat).integrate('control_region',s).dim() != 2: continue
            plot.plot1d(hists[key].integrate('category',cat).integrate('control_region',s), overlay="process", ax=ax, clear=False, stack=True, fill_opts=fill_opts, error_opts=error_opts)
            if str(s) != "signal":
                plot.plot1d(data_hists[key].integrate('category',cat).integrate('control_region',s), overlay="process", ax=ax, clear=False, error_opts=data_err_opts)
                plot.plotratio(data_hists[key].integrate('category',cat).integrate('control_region',s).integrate('process'), hists[key].integrate('category',cat).integrate('control_region',s).sum('process'), ax=rax, overflow='over', error_opts=data_err_opts, denom_fill_opts={}, guide_opts={}, unc='num')
            rax.set_ylabel('Ratio (data/MC)')
            rax.set_ylim(0.5,1.5)
            ax._get_lines.prop_cycler = ax._get_patches_for_fill.prop_cycler
            args = {'linestyle':'--', 'linewidth': 5}
            height = 2*max(max(hists[key].integrate('category',cat).integrate('control_region',s).integrate('process').values()[()]),max(data_hists[key].integrate('category',cat).integrate('control_region',s).integrate('process').values()[()]))
            ax.autoscale(axis='x', tight=True)
            ax.set_yscale('log')
            ax.set_ylim(bottom=0.1, top=height)
            leg = ax.legend()
            coffee = plt.text(0., 1., u"☕ "+cat_names[str(cat)]+" "+str(s), fontsize=20, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
            lumi = plt.text(1., 1., r"59.97 fb$^{-1}$ (13 TeV)", fontsize=20, horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)
            plot_path = os.path.abspath('stack') 
            plot_name = key+'_'+str(cat)+'_'+str(s)+'_stack.png'
            fig.savefig(os.path.join(plot_path, plot_name))
    
if not os.path.exists('unstack'): os.makedirs('unstack')
for key in hists.keys():
    if key=='sumw': continue
    
    for cat in hists[key].identifiers('category'):
        for s in hists[key].integrate('category',cat).identifiers('control_region'):
            if str(s)=="QCD": continue
            #pdb.set_trace()
            args = {'linestyle':'--','linewidth':2}
            fig, ax = plt.subplots(1, 1, figsize=(10,10))
            ax.set_prop_cycle(cycler(color=colors))
            if hists[key].integrate('category',cat).integrate('control_region',s).dim() != 2: continue
            plot.plot1d(hists[key].integrate('category',cat).integrate('control_region',s), ax=ax, overlay="process", clear=False, stack=False, line_opts={}, density=1)
            ax.autoscale(axis='x', tight=True)
            leg = ax.legend()
            coffee = plt.text(0., 1., u"☕ "+cat_names[str(cat)]+" "+str(s), fontsize=28, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
            lumi = plt.text(1., 1., r"1 fb$^{-1}$ (13 TeV)", fontsize=16, horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)
            plot_path = os.path.abspath('unstack')
            plot_name = key+'_'+str(cat)+'_'+str(s)+'.png'
            fig.savefig(os.path.join(plot_path, plot_name))
            
# for key in hists.keys():
#     if key=='sumw': continue

#     fig, ax = plt.subplots(1, 1, figsize=(7,7))
#     ax.set_prop_cycle(cycler(color=colors))
#     plot.plot1d(hists[key].project('control_region').project('category'), ax=ax, overlay="process", clear=False, stack=False, line_opts={},density=1)
#     ax.autoscale(axis='x', tight=True)
#     ax.set_yscale('log')
#     ax.set_ylim(.1, None)
#     leg = ax.legend()
#     coffee = plt.text(0., 1., u"☕ Test", fontsize=28, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
#     #lumi = plt.text(1., 1., r"1 fb$^{-1}$ (13 TeV)", fontsize=16, horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)

# fig, ax = plt.subplots(1, 1, figsize=(7,7))
# ax.set_prop_cycle(cycler(color=colors))
# plot.plot1d(hists["recoil"].project('control_region').project('category'), ax=ax, overlay="process", clear=False, stack=False, line_opts={},density=1)
# ax.autoscale(axis='x', tight=True)
# ax.set_yscale('log')
# leg = ax.legend()
# coffee = plt.text(0., 1., u"☕", fontsize=28, horizontalalignment='left', verticalalignment='bottom', transform=ax.transAxes)
# lumi = plt.text(1., 1., r"1 fb$^{-1}$ (130 TeV)", fontsize=16, horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes)
