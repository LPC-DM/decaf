import uproot as up
import numpy as np
import awkward as ak
import matplotlib.pyplot as plt
import mplhep as hep
from optparse import OptionParser
from cycler import cycler
import itertools

hep.style.use("CMS")

#### Define region, year, signal, and recoil bin ####
regions = {
        "signal": "sr",
        "singleew": "wecr",
        "singlemw": "wmcr",
        "singlemtop":"tmcr",
        "singleetop":"tecr",
}

years = {
        "2016": "2016",
        "2017": "2017",
        "2018": "2018"
}

lumi = {
        "2016": 35.92,
        "2017": 41.53,
        "2018": 59.74
}

category = {
        "pass": "pass",
        "fail": "fail"
}

recoilDict = {
        "recoil0": "recoil0",
        "recoil1": "recoil1",
        "recoil2": "recoil2",
        "recoil3": "recoil3",
        "recoil4": "recoil4"
}

processNames = {
    'qcdMC': 'QCD',
    'tt': r'$t\bar{t}$',
    'ttMC': r'$t\bar{t}$',
    'stMC': 'Single t',
    'vvMC': 'Diboson',
    'hbbMC': r'$H\rightarrow b\bar{b}$',
    'dyjetsMC': 'DY+jets',
    'wjets': 'W+jets',
    'wjetsMC': 'W+jets',
    'zjets': 'Z+jets'
    #'Mhs_50': 'Signal'
}

colorDict = {
    'DY+jets': '#666666',
    r'$H\rightarrow b\bar{b}$': '#a6761d',
    'QCD': '#e6ab02',
    'Diboson': '#66a61e',
    'Single t': '#e7298a',
    r'$t\bar{t}$': '#7570b3',
    'W+jets': '#d95f02',
    'Z+jets':'#1b9e77'
}

#extralabels = {
#        "singlemw":"Single-#mu b-vetoed CR",
#        "singlemtop":"Single-#mu b-tagged CR",
#        "signal":"Signal region",
#        "singleetop":"Single-e b-tagged CR",
#        "singleew":"Single-e b-vetoed CR"
#}

def plotTCR(region, year, recoil):
    ### Make printing which region, year, signal, and recoil are considered ###
    print('-----------------------------------------------------------')
    print('You are considering region:', regions[region])
    print('Considering year is:', years[year])
    print('Pass category')
    print('Recoil bin is:', recoilDict[recoil])
    print('Opening', regions[region]+years[year]+"pass"+recoilDict[recoil])
    print('-----------------------------------------------------------')

    ### Define processes per region ###
    processesT = ['dyjetsMC', 'hbbMC', 'qcdMC', 'vvMC', 'stMC', 'tt', 'wjetsMC']
    t_label = ['DY+jets', r'$H\rightarrow b\bar{b}$', 'QCD', 'Diboson', 'Single t', r'$t\bar{t}$', 'W+jets']

    processes = processesT
    mc_labels = t_label

    ### Check directory ###
    f = up.open("postfitshapes.result.Mz1000mhs70Mdm150.mass120to300.root")

    prefit_dir2 = regions[region]+years[year]+"passmass120to300"+recoilDict[recoil]+"_prefit"
    postfit_dir2 = regions[region]+years[year]+"passmass120to300"+recoilDict[recoil]+"_postfit"

    ### First, using prefit histogram to define histogram edges
    dir2 = f[prefit_dir2]
    prefit_bins, edges2 = dir2["TotalProcs"].to_numpy()

    ### To add comma in the array, convert list and back to numpy
    prefit = prefit_bins.tolist()
    edges = edges2.tolist()

    fig, (ax, rax) = plt.subplots(2, 1, figsize=(10,10), gridspec_kw=dict(height_ratios=[3, 1], hspace=0.07), sharex=True)
    errps = {'hatch':'////', 'facecolor':'none', 'lw': 0, 'color': 'k', 'alpha': 0.4}
    ax.set_ylabel('Events/GeV')
    ax._get_lines.prop_cycler = ax._get_patches_for_fill.prop_cycler
    args = {'linestyle':'--', 'linewidth': 5}
    ax.set_yscale('log')
    ax.set_ylim(1e-2, 5e+5)

    if year == "2016":
        hep.cms.label(ax=ax, loc=0, lumi=36, year=2016)
    elif year == "2017":
        hep.cms.label(ax=ax, loc=0, lumi=42, year=2017)
    elif year == "2018":
        hep.cms.label(ax=ax, loc=0, lumi=60, year=2018)

    ### Move to draw postfit
    sum_postfit = np.zeros(prefit_bins.size)
    dir2 = f[postfit_dir2]

    process_bin = []
    keys = [x.split(';')[0] for x in dir2.keys()]

    for j in processes:
        print('Which process', j)
        if not j in keys:
            print('Not found, skip this process')
            print('Remove its label too \n')
            mc_labels.remove(processNames[j])
            continue
        bins2, _ = dir2[j].to_numpy()
        sum_postfit += bins2
        process_bin.append(bins2.tolist())

    postfit = sum_postfit.tolist()
    colors = [colorDict[x] for x in mc_labels]
    ax.set_prop_cycle(cycler(color=colors))

    ### Try to draw stack plots
    hep.histplot(process_bin, edges, ax=ax, stack=True, histtype='fill', edgecolor = 'k', linewidth=1, label=mc_labels)

    ### Draw Stat. unc.
    ax.stairs(
        values=sum_postfit + np.sqrt(sum_postfit),
        baseline=sum_postfit - np.sqrt(sum_postfit),
        edges=edges, **errps, label='Stat. unc.')

    hep.histplot(postfit, edges, ax=ax, label=["SM total (post-fit)"], color='b', linewidth=3)
    hep.histplot(prefit, edges, ax=ax, label=["SM total (pre-fit)"], color='r', linestyle='dashed', linewidth=2)

    ### Call data ###
    data2, _ = dir2["data_obs"].to_numpy()

    hep.histplot(data2, edges, ax=ax, histtype='errorbar', label="Data", color='k')

    from hist.intervals import ratio_uncertainty
    yerr = ratio_uncertainty(data2, sum_postfit, 'poisson')
    rax.stairs(1+yerr[1], edges=edges, baseline=1-yerr[0], **errps)

    hep.histplot(data2/prefit_bins, edges, yerr=np.sqrt(data2)/prefit_bins, ax=rax, histtype='errorbar', color='r', capsize=4, label="Prefit")
    hep.histplot(data2/sum_postfit, edges, yerr=np.sqrt(data2)/sum_postfit, ax=rax, histtype='errorbar', color='b', capsize=4, label="Postfit")

    rax.axhline(1, ls='--', color='k')
    rax.set_ylim(0.5, 1.5)
    rax.set_xlabel('$p_{T}^{miss}$ [GeV]')
    rax.set_ylabel('Obs/Exp', fontsize=15, loc='center')
    rax.legend(loc='upper right', fontsize=12, ncol=2)

    name = regions[region]+years[year]+"pass"+recoilDict[recoil]
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc='upper right', fontsize=12, ncol=2)
    fig.savefig('plots/'+name+'.png')
    plt.close('all')

def plotWCR(region, year, cate, recoil):

    #extralabel = extralabels[region]+' '+category[cate]+' '+recoilDict[recoil]

    ### Make printing which region, year, signal, and recoil are considered ###
    print('-----------------------------------------------------------')
    print('You are considering region:', regions[region])
    print('Considering year is:', years[year])
    print('Signal process is:', category[cate])
    print('Recoil bin is:', recoilDict[recoil])
    print('Opening', regions[region]+years[year]+category[cate]+recoilDict[recoil])
    print('-----------------------------------------------------------')

    ### Define processes per region ###
    processesW1 = ['dyjetsMC', 'hbbMC', 'qcdMC', 'vvMC', 'stMC', 'ttMC', 'wjets']
    processesW2 = ['dyjetsMC', 'hbbMC', 'qcdMC', 'vvMC', 'stMC', 'tt', 'wjets']
    w_label = ['DY+jets', r'$H\rightarrow b\bar{b}$', 'QCD', 'Diboson', 'Single t', r'$t\bar{t}$', 'W+jets']

    if category[cate] == "pass":
        processes1 = processesW1
        processes2 = processesW2
        mc_labels = w_label

    if category[cate] == "fail":
        processes1 = processesW1
        processes2 = processesW1
        mc_labels = w_label

    print('List of processes', processes1, processes2, '\n')

    ### Check directory ###
    f1 = up.open("postfitshapes.result.Mz1000mhs70Mdm150.mass40to120.root")
    f2 = up.open("postfitshapes.result.Mz1000mhs70Mdm150.mass120to300.root")

    prefit_dir1 = regions[region]+years[year]+category[cate]+"mass40to120"+recoilDict[recoil]+"_prefit"
    prefit_dir2 = regions[region]+years[year]+category[cate]+"mass120to300"+recoilDict[recoil]+"_prefit"

    postfit_dir1 = regions[region]+years[year]+category[cate]+"mass40to120"+recoilDict[recoil]+"_postfit"
    postfit_dir2 = regions[region]+years[year]+category[cate]+"mass120to300"+recoilDict[recoil]+"_postfit"

    ### First, using prefit histogram to define histogram edges
    dir1 = f1[prefit_dir1]
    dir2 = f2[prefit_dir2]

    bins1, edges1 = dir1["TotalProcs"].to_numpy()
    bins2, edges2 = dir2["TotalProcs"].to_numpy()
    edges2 = np.delete(edges2, 0) ## remove duplicate number

    prefit_bin = np.concatenate((bins1, bins2), axis=None)
    total_edges = np.concatenate((edges1, edges2), axis=None)

    ### To add comma in the array, convert list and back to numpy
    prefit = prefit_bin.tolist()
    edges = total_edges.tolist()

    #fig, ax = plt.subplots(figsize=(8, 8), constrained_layout=True)
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(10,10), gridspec_kw=dict(height_ratios=[3, 1], hspace=0.07), sharex=True)
    errps = {'hatch':'////', 'facecolor':'none', 'lw': 0, 'color': 'k', 'alpha': 0.4}
    ax.set_ylabel('Events/GeV')
    ax._get_lines.prop_cycler = ax._get_patches_for_fill.prop_cycler
    args = {'linestyle':'--', 'linewidth': 5}
    ax.set_yscale('log')
    ax.set_ylim(1e-2, 5e+5)

    if year == "2016":
        hep.cms.label(ax=ax, loc=0, lumi=36, year=2016)
    elif year == "2017":
        hep.cms.label(ax=ax, loc=0, lumi=42, year=2017)
    elif year == "2018":
        hep.cms.label(ax=ax, loc=0, lumi=60, year=2018)

    ### Move to draw postfit
    sum_postfit = np.zeros(prefit_bin.size)

    dir1 = f1[postfit_dir1]
    dir2 = f2[postfit_dir2]

    process_bin = []

    keys1 = [x.split(';')[0] for x in dir1.keys()]
    keys2 = [x.split(';')[0] for x in dir2.keys()]

    #for process in processes:
    for (i, j) in zip(processes1, processes2):
        print('Which process', i, j)
        if not i in keys1:
            tbins1 = np.zeros(bins1.size)
        else:
            tbins1, _ = dir1[i].to_numpy()

        if not j in keys2:
            tbins2 = np.zeros(bins2.size)
        else:
            tbins2, _ = dir2[j].to_numpy()

        postfit_bin = np.concatenate((tbins1, tbins2), axis=None)
        if postfit_bin.sum() == 0.:
            print('Not found, skip this process')
            print('Remove its label too \n')
            mc_labels.remove(processNames[i])
            continue
        sum_postfit += postfit_bin
        process_bin.append(postfit_bin.tolist())

    postfit = sum_postfit.tolist()
    colors = [colorDict[x] for x in mc_labels]
    ax.set_prop_cycle(cycler(color=colors))

    ### Try to draw stack plots
    hep.histplot(process_bin, edges, ax=ax, stack=True, histtype='fill', edgecolor = 'k', linewidth=1, label=mc_labels)

    ### Draw Stat. unc.
    ax.stairs(
        values=sum_postfit + np.sqrt(sum_postfit),
        baseline=sum_postfit - np.sqrt(sum_postfit),
        edges=edges, **errps, label='Stat. unc.')

    hep.histplot(postfit, edges, ax=ax, label=["SM total (post-fit)"], color='b', linewidth=3)
    hep.histplot(prefit, edges, ax=ax, label=["SM total (pre-fit)"], color='r', linestyle='dashed', linewidth=2)

    ### Call data ###
    data1, _ = dir1["data_obs"].to_numpy()
    data2, _ = dir2["data_obs"].to_numpy()
    data = np.concatenate((data1, data2), axis=None)

    hep.histplot(data, edges, ax=ax, histtype='errorbar', label="Data", color='k')

    from hist.intervals import ratio_uncertainty
    yerr = ratio_uncertainty(data, sum_postfit, 'poisson')
    rax.stairs(1+yerr[1], edges=edges, baseline=1-yerr[0], **errps)

    hep.histplot(data/prefit_bin, edges, yerr=np.sqrt(data)/prefit_bin, ax=rax, histtype='errorbar', color='r', capsize=4, label="Prefit")
    hep.histplot(data/sum_postfit, edges, yerr=np.sqrt(data)/sum_postfit, ax=rax, histtype='errorbar', color='b', capsize=4, label="Postfit")

    rax.axhline(1, ls='--', color='k')
    rax.set_ylim(0.5, 1.5)
    rax.set_xlabel('$p_{T}^{miss}$ [GeV]')
    rax.set_ylabel('Obs/Exp', fontsize=15, loc='center')
    rax.legend(loc='upper right', fontsize=12, ncol=2)

    name = regions[region]+years[year]+category[cate]+recoilDict[recoil]
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc='upper right', fontsize=12, ncol=2)
    fig.savefig('plots/'+name+'.png')
    plt.close('all')

def plotSR(region, year, cate, recoil):

    #extralabel = extralabels[region]+' '+category[cate]+' '+recoilDict[recoil]

    ### Make printing which region, year, signal, and recoil are considered ###
    print('-----------------------------------------------------------')
    print('You are considering region:', regions[region])
    print('Considering year is:', years[year])
    print('Signal process is:', category[cate])
    print('Recoil bin is:', recoilDict[recoil])
    print('Opening', regions[region]+years[year]+category[cate]+recoilDict[recoil])
    print('-----------------------------------------------------------')

    ### Define processes per region ###
    processesS1 = ['dyjetsMC', 'hbbMC', 'qcdMC', 'vvMC', 'stMC', 'ttMC', 'wjets', 'zjets']
    processesS2 = ['dyjetsMC', 'hbbMC', 'qcdMC', 'vvMC', 'stMC', 'tt', 'wjets', 'zjets']
    s_label = ['DY+jets', r'$H\rightarrow b\bar{b}$', 'QCD', 'Diboson', 'Single t', r'$t\bar{t}$', 'W+jets', 'Z+jets']

    processesS_r4 = ['dyjetsMC', 'hbbMC', 'qcdMC', 'vvMC', 'stMC', 'ttMC', 'wjetsMC', 'zjets']

    if category[cate] == "pass":
        processes1 = processesS1
        processes2 = processesS2
        if recoilDict[recoil] == "recoil4":
            processes1 = processesS_r4
            processes2 = processesS_r4
        mc_labels = s_label

    if category[cate] == "fail":
        processes1 = processesS1
        processes2 = processesS1
        mc_labels = s_label

    print('List of processes', processes1, processes2, '\n')

    ### Check directory ###
    f1 = up.open("postfitshapes.result.Mz1000mhs70Mdm150.mass40to120.root")
    f2 = up.open("postfitshapes.result.Mz1000mhs70Mdm150.mass120to300.root")

    prefit_dir1 = regions[region]+years[year]+category[cate]+"mass40to120"+recoilDict[recoil]+"_prefit"
    prefit_dir2 = regions[region]+years[year]+category[cate]+"mass120to300"+recoilDict[recoil]+"_prefit"

    postfit_dir1 = regions[region]+years[year]+category[cate]+"mass40to120"+recoilDict[recoil]+"_postfit"
    postfit_dir2 = regions[region]+years[year]+category[cate]+"mass120to300"+recoilDict[recoil]+"_postfit"

    ### First, using prefit histogram to define histogram edges
    dir1 = f1[prefit_dir1]
    dir2 = f2[prefit_dir2]

    bins1, edges1 = dir1["TotalProcs"].to_numpy()
    bins2, edges2 = dir2["TotalProcs"].to_numpy()
    edges2 = np.delete(edges2, 0) ## remove duplicate number

    prefit_bin = np.concatenate((bins1, bins2), axis=None)
    total_edges = np.concatenate((edges1, edges2), axis=None)

    ### To add comma in the array, convert list and back to numpy
    prefit = prefit_bin.tolist()
    edges = total_edges.tolist()

    #fig, ax = plt.subplots(figsize=(8, 8), constrained_layout=True)
    fig, (ax, rax) = plt.subplots(2, 1, figsize=(10,10), gridspec_kw=dict(height_ratios=[3, 1], hspace=0.07), sharex=True)
    errps = {'hatch':'////', 'facecolor':'none', 'lw': 0, 'color': 'k', 'alpha': 0.4}
    ax.set_ylabel('Events/GeV')
    ax._get_lines.prop_cycler = ax._get_patches_for_fill.prop_cycler
    args = {'linestyle':'--', 'linewidth': 5}
    ax.set_yscale('log')
    ax.set_ylim(1e-2, 5e+5)

    if year == "2016":
        hep.cms.label(ax=ax, loc=0, lumi=36, year=2016)
    elif year == "2017":
        hep.cms.label(ax=ax, loc=0, lumi=42, year=2017)
    elif year == "2018":
        hep.cms.label(ax=ax, loc=0, lumi=60, year=2018)

    ### Move to draw postfit
    sum_postfit = np.zeros(prefit_bin.size)

    dir1 = f1[postfit_dir1]
    dir2 = f2[postfit_dir2]

    process_bin = []

    keys1 = [x.split(';')[0] for x in dir1.keys()]
    keys2 = [x.split(';')[0] for x in dir2.keys()]

    #for process in processes:
    for (i, j) in zip(processes1, processes2):
        print('Which process', i, j)
        if not i in keys1:
            tbins1 = np.zeros(bins1.size)
        else:
            tbins1, _ = dir1[i].to_numpy()

        if not j in keys2:
            tbins2 = np.zeros(bins2.size)
        else:
            tbins2, _ = dir2[j].to_numpy()
        postfit_bin = np.concatenate((tbins1, tbins2), axis=None)

        if postfit_bin.sum() == 0.:
            print('Not found, skip this process')
            print('Remove its label too \n')
            mc_labels.remove(processNames[i])
            continue

        sum_postfit += postfit_bin
        process_bin.append(postfit_bin.tolist())

    postfit = sum_postfit.tolist()
    colors = [colorDict[x] for x in mc_labels]
    ax.set_prop_cycle(cycler(color=colors))

    ### Try to draw stack plots
    hep.histplot(process_bin, edges, ax=ax, stack=True, histtype='fill', edgecolor = 'k', linewidth=1, label=mc_labels)

    ### Draw Stat. unc.
    ax.stairs(
        values=sum_postfit + np.sqrt(sum_postfit),
        baseline=sum_postfit - np.sqrt(sum_postfit),
        edges=edges, **errps, label='Stat. unc.')

    hep.histplot(postfit, edges, ax=ax, label=["SM total (post-fit)"], color='b', linewidth=3)
    hep.histplot(prefit, edges, ax=ax, label=["SM total (pre-fit)"], color='r', linestyle='dashed', linewidth=2)

    ### Call data ###
    data1, _ = dir1["data_obs"].to_numpy()
    data2, _ = dir2["data_obs"].to_numpy()
    data = np.concatenate((data1, data2), axis=None)

    hep.histplot(data, edges, ax=ax, histtype='errorbar', label="Data", color='k')

    from hist.intervals import ratio_uncertainty
    yerr = ratio_uncertainty(data, sum_postfit, 'poisson')
    rax.stairs(1+yerr[1], edges=edges, baseline=1-yerr[0], **errps)

    hep.histplot(data/prefit_bin, edges, yerr=np.sqrt(data)/prefit_bin, ax=rax, histtype='errorbar', color='r', capsize=4, label="Prefit")
    hep.histplot(data/sum_postfit, edges, yerr=np.sqrt(data)/sum_postfit, ax=rax, histtype='errorbar', color='b', capsize=4, label="Postfit")

    rax.axhline(1, ls='--', color='k')
    rax.set_ylim(0.5, 1.5)
    rax.set_xlabel('$p_{T}^{miss}$ [GeV]')
    rax.set_ylabel('Obs/Exp', fontsize=15, loc='center')
    rax.legend(loc='upper right', fontsize=12, ncol=2)

    name = regions[region]+years[year]+category[cate]+recoilDict[recoil]
    handles, labels = ax.get_legend_handles_labels()
    ax.legend(handles[::-1], labels[::-1], loc='upper right', fontsize=12, ncol=2)
    fig.savefig('plots/'+name+'.png')
    plt.close('all')

#### Call Functions ####
if __name__ == "__main__":

    bins = ['recoil0', 'recoil1', 'recoil2', 'recoil3', 'recoil4']

    parser = OptionParser()
    parser.add_option("-y", "--year", help="year", dest="year", default="")
    (options, args) = parser.parse_args()
    year = options.year

    for iregion in ['singlemw', 'singleew']:
        for ica in ["pass", "fail"]:
            for ibin in bins:
                #plotWCR(iregion, year, ica, ibin)
                try:
                    plotWCR(iregion, year, ica, ibin)
                except:
                    print("Directory do not exist in %s file! \n" % ("yes"))
                    pass

    for iregion in ['signal']:
        for ica in ["pass", "fail"]:
            for ibin in bins:
                plotSR(iregion, year, ica, ibin)
                #try:
                #    plotSR(iregion, year, ica, ibin)
                #except:
                #    print("Directory do not exist in %s file! \n" % ("yes"))
                #    pass

    for iregion in ['singlemtop', 'singleetop']:
        for ibin in ['recoil0', 'recoil1', 'recoil2', 'recoil3']:
            plotTCR(iregion, year, ibin)
            #try:
            #    plotTCR(iregion, year, ibin)
            #except:
            #    print("Directory do not exist in %s file! \n" % ("yes"))
            #    pass
