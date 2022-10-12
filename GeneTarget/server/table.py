#Generate the table.csv = dataframe
from os import listdir
from os.path import isfile, join
mypath='events/'

data={}
genes=set()
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]
for fileName in onlyfiles:
  i=fileName.find("_")
gene=fileName[:i]
sample=fileName[i+1:]
genes.add(gene)
lines= [l.strip() for l in open('events/'+fileName)]
svs=filter(lambda l: len(l)>0 and l[0]!='#',lines)
numEvents=len(list(svs))
if sample not in data:
  data[sample]={}
data[sample][gene]=numEvents

genes=list(genes)
outputFile=open("out.csv",'w')
header=[sample] + genes
outputFile.write(",".join(header)+"\n")

for sample, events in data.items():
  line=[sample]+[str(events[g]) for g in genes]
outputFile.write(",".join(line)+"\n")
