#### Make text2workspace command ####
signals = ['Mz2500_mhs90_Mdm1250', 'Mz500_mhs70_Mdm150', 'Mz2000_mhs50_Mdm1000', 'Mz200_mhs50_Mdm150', 'Mz2000_mhs50_Mdm500', 'Mz2500_mhs50_Mdm750', 'Mz3000_mhs70_Mdm1000', 'Mz300_mhs70_Mdm100', 'Mz500_mhs90_Mdm150', 'Mz200_mhs90_Mdm100', 'Mz1000_mhs50_Mdm500', 'Mz3000_mhs90_Mdm1500', 'Mz200_mhs70_Mdm100', 'Mz2500_mhs90_Mdm750', 'Mz300_mhs50_Mdm150', 'Mz2500_mhs50_Mdm1250', 'Mz1000_mhs50_Mdm1000', 'Mz500_mhs50_Mdm150', 'Mz500_mhs90_Mdm500', 'Mz200_mhs70_Mdm150', 'Mz300_mhs50_Mdm100', 'Mz2000_mhs90_Mdm1500', 'Mz500_mhs70_Mdm500', 'Mz1000_mhs90_Mdm1000', 'Mz2000_mhs70_Mdm1500', 'Mz3000_mhs70_Mdm1500', 'Mz3000_mhs50_Mdm1500', 'Mz500_mhs50_Mdm250', 'Mz300_mhs90_Mdm150', 'Mz300_mhs70_Mdm150', 'Mz1000_mhs70_Mdm150', 'Mz2500_mhs70_Mdm1250', 'Mz3000_mhs50_Mdm1000', 'Mz2000_mhs90_Mdm500', 'Mz200_mhs50_Mdm100', 'Mz3000_mhs90_Mdm1000', 'Mz500_mhs90_Mdm250', 'Mz2000_mhs50_Mdm1500', 'Mz300_mhs90_Mdm100', 'Mz1000_mhs90_Mdm500', 'Mz1000_mhs90_Mdm150', 'Mz2000_mhs70_Mdm1000', 'Mz1000_mhs70_Mdm500', 'Mz500_mhs70_Mdm250', 'Mz2000_mhs70_Mdm500', 'Mz1000_mhs50_Mdm150', 'Mz2500_mhs70_Mdm750', 'Mz500_mhs50_Mdm500', 'Mz1000_mhs70_Mdm1000', 'Mz200_mhs90_Mdm150', 'Mz2000_mhs90_Mdm1000']

tmp = []
for signal in signals:
    label = ' --PO map=.*/'+signal+':r_'+signal+'[1,0,2]'
    tmp.append(label)
option = ''.join(tmp)

for year in ['2016', '2017', '2018']:
    command = "text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel --PO verbose"+option+' Darkhiggs'+year+'.txt -o Darkhiggs'+year+'.root'
    print(command, '\n')
