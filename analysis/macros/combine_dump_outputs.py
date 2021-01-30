import ROOT as rt
import glob, os
import sys

##### example: plots/darkhiggs2018/dump_postfit ####
dirname = sys.argv[1]
print(dirname.split('/')[2].split('_')[1])
flag = dirname.split('/')[2].split('_')[1]

def makeinputs(wsname):
    filelist = glob.glob(dirname+'/*'+wsname+'*.root')
    print(filelist)
    return filelist

def histcombine(region, year, signalflag, recoil):
    #### Define region, year, signal, and recoil bin ####
    darkhiggs_regions = {
            "signal": "sr",
            "singleetop": "tecr",
            "singlemtop": "tmcr",
            "singleew": "wecr",
            "singlemw": "wmcr",
            "photon": "gcr",
            "diele": "zecr",
            "dimu": "zmcr"
    }

    years = {
            "2016": "2016",
            "2017": "2017",
            "2018": "2018"
    }

    signalprocess = {
            "mhs": "pass",
            "mjet": "fail"
    }

    recoilbin = {
            "recoil0": "recoil0",
            "recoil1": "recoil1",
            "recoil2": "recoil2",
            "recoil3": "recoil3",
            "recoil4": "recoil4"
    }

    workspace = darkhiggs_regions[region]+years[year]+signalprocess[signalflag]+recoilbin[recoil]
    inputlist = makeinputs(workspace)

    if not inputlist:
        raise RuntimeError

    outputdir = 'plots/darkhiggs2016/dump'
    if not os.path.exists(outputdir):
        os.makedirs(outputdir)

    foutput = rt.TFile(outputdir+'/'+flag+'_'+workspace+'.root', 'recreate')
    foutput.mkdir(workspace+'_'+flag)

    for ifile in inputlist:
        process = ''
        if 'Mhs_50_morph' in str(ifile):
            process = 'Mhs_50'
        elif 'data' in str(ifile):
            process = 'data'
        elif 'morph' in str(ifile):
            process = ifile.split('_')[-2]
        else:
            process = ifile.split('/')[3].split('_')[1]

        fin = rt.TFile.Open(ifile)
        if fin == None:
            print('Cannot open the file!')

        if process != 'data':
            hist = fin.Get('hist_'+ifile.split('/')[-1].split('.')[0]+'__fjmass')
        else:
            print('hist_shapeBkg__data_'+workspace+'Pdf__fjmass')
            hist = fin.Get('hist_shapeBkg__data_'+workspace+'Pdf__fjmass')
        foutput.cd(workspace+'_'+flag)
        hist.Write(process)
        fin.Close()

    foutput.Close()

if __name__ == '__main__':
    dh_regions = ['singlemw', 'singlemtop', 'singleew', 'singleetop', 'photon', 'diele', 'dimu', 'signal']
    sigs = ['mhs', 'mjet']
    bins = ['recoil0', 'recoil1', 'recoil2', 'recoil3', 'recoil4']

    for iregion in dh_regions:
        for isig in sigs:
            for ibin in bins:
                try:
                    histcombine(iregion, '2016', isig, ibin)
                except:
                    print('There are no inputs for %s %s %s \n' %(iregion, isig, ibin))
                    pass
