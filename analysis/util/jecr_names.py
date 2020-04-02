import os


def change_name(directory,types=None):
    for filename in os.listdir(directory):
        if '~' in filename: continue
        if 'junc' in filename: continue
        if 'jec' in filename: continue
        if 'jr' in filename: continue
        if 'jersf' in filename: continue
        new_filename=''
        if 'Uncertainty' in filename: new_filename = filename.split('.')[0]+'.junc.txt'
        if 'SF' in filename: new_filename = filename.split('.')[0]+'.jersf.txt'
        if 'jr' in directory: new_filename = filename.split('.')[0]+'.jr.txt' 
        if 'jec' in directory: new_filename = filename.split('.')[0]+'.jec.txt'
        if types is not None: new_filename = filename.split('.')[0]+'.'+types+'.txt'
        #print('mv '+directory+'/'+filename+' '+directory+'/'+new_filename)
        os.system('mv '+directory+'/'+filename+' '+directory+'/'+new_filename)

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-f', '--folder', help='folder', dest='folder')
    parser.add_option('-t', '--types', help='types', dest='types')
    (options, args) = parser.parse_args()
    change_name(options.folder)
    if options.types: change_name(options.folder,options.types) 
