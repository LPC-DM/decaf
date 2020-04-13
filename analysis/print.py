from coffea.util import load

hists = load('hists/monotop2018.merged')
#print(hists['recoil'].identifiers('dataset'))

dataset = []
for p in hists['recoil'].identifiers('dataset'):
    dataset.append(str(p))

print(*dataset, sep = "\n")
