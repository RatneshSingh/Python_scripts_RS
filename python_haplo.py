#!/usr/bin/python
import sys
#print "More than"+sys.argv[2]+" missing info"
#print sys.argv[1]
try:
  fhandle=open(sys.argv[1],"r")
except:
  print "unable to open",sys.argv[1]
  exit
print "Catalog_Id\tCnt\t",sys.argv[1];
headers=[]
alpersamp={}
for haplo in fhandle:
  alleles={}
  if haplo.startswith("Cat"):
    headers=haplo.split("\t")
    for hkey in headers:alpersamp={hkey:{0:0}}
    continue
  #if "-" in haplo:continue
  #print "This line has",haplo.count('-')," missing data"+haplo
  if haplo.count('-') > int(sys.argv[2]):
    #print "More than"+sys.argv[2]+"missing info"
    continue
  #print "Reading line",haplo
  line=haplo.split()
  #print("Lenght of line",len(line))
  for i in range(2,len(line)-1):
    #if "consensus" in line[i]:continue
    if "-" not in line[i]:
      if alpersamp[headers[i]][line[i].count("/") + 1 ] < 1:alpersamp={headers[i]:{line[i].count("/")+1:1}}
      else: alpersamp[headers[i]][line[i].count("/") + 1 ]=alpersamp[headers[i]][line[i].count("/") + 1 ]+1
      for j in line[i].split("/"):
        alleles[j]=1

  big_allele="/".join(alleles.keys())
  if len(big_allele) < 1:big_allele="-"
  #print "printing after join",big_allele
  #print "\t".join([line[0],line[1],big_allele,str(len(alleles))])


## print allele num in each sample
for samp in alpersamp:
  for allnum in samp:
    print "\t".join([samp,allnum,alpersamp[samp][allnum]])
