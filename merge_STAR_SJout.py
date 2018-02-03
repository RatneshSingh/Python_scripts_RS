#!/usr/bin/python
from sys import argv, exit
from BCBio import GFF
#from collections import defaultdict
# from BCBio.GFF import GFFExaminer
import pprint
import getopt

### functions

def usage(mesg=" "):
    
    print("\n*****  {}  ****\n".format(mesg))
    help="""
    python3 {} -f SJ.out.tab1,SJ.out.tab2,SJ.out.tab....N  -g GFF3 formatted annotation file
    options:
    -o	Out file to save results
    -c  type of intron motif to use [0]
    	0: All types
    	1: GT/AG
    	2: CT/AC
    	3: GC/AG
    	4: CT/GC
    	5: AT/AC
    	6: GT/AT 
    -m  Report all (uniquely +  multi-mapped) reads across junction [False]
print("\n*****  {}  ****\n".format(mesg))
        By default script reports uniquely mapped reads only.
    """.format(argv[0],)
    print(help)
    exit()

def overlaps(frst,scnd):
    frst=sorted(frst)
    scnd=sorted(scnd)
    return min(frst[-1],scnd[-1]) - max(frst[0],scnd[0]) > 0



spmot={0:'All',1:'GT/AG', 2: 'CT/AC',3: 'GC/AG',4: 'CT/GC',5: 'AT/AC',6: 'GT/AT'}



try:
    opts, args = getopt.getopt(argv[1:], 'hmc:g:f:o:')
except:
    print("One of the argument in:", argv, " is invalid:")
    usage()

files = str()
canon=0
multi=0
for opt, arg in opts:
    if opt == "-h":
        usage()
    elif opt in ("-g"):
        gff_file = arg
        print("Found gff file:", gff_file)
    elif opt in ("-f"):
        files = arg.split(',')
        #print("Found files:", files)
    elif opt in ("-o"):
        ofile = arg
        print("Results will be saved in", ofile)
    elif opt in ("-c"):
        canon = int(arg)
        print("Results will include introns with ", spmot[canon]," motifs")
    elif opt in ("-m"):
        multi = 1
        print("Results will include multimapped reads in read count")
    else:
        print("Use legal flags\n", argv)
        print("error:", args)
        usage("Unable to parse argument")
        

## 
try:gff_file
except: usage("No gff file was provided")

if len(files) > 0: print("\nTotal number of SJ out files ".format(len(files)))
else: usage("No SJ.out.tab file was provided")




sj={}
gstrt={}



## read file to collect sj boundry and read count
for file in files:
    sj[file]={}
    print("Reading ",file)
    f=open(file,'r')
    for line in f:
        elem=[]
        elem=line.split()

        if canon > 0 and int(elem[4]) != canon: continue
        ### calculate mapped read uniq or multi + unique
        if multi < 1:
            mread=int(elem[6])
        else:
            mread=int(elem[6])+int(elem[7])

        try:
            sj[file][elem[0]].update({int(elem[1]):{int(elem[2]):{'uread':mread}}})
        except:
            sj[file].update({elem[0]:{int(elem[1]):{int(elem[2]):{'uread':mread}}}})

        try:
            gstrt[elem[0]].update({int(elem[1]):{int(elem[2]):1}})
        except:
            gstrt[elem[0]]={int(elem[1]):{int(elem[2]):1}}

    f.close()

#print sj
print("Number of files read:",str(len(sj)))
#pprint.pprint(sj)
#quit()

#### collect gene name information from gff file to identify exons later.
g=open(gff_file,"r")
gene_info={}
for gline in g:
    print("Processing Line" + gline)
    if gline.startswith("#"): continue
    if gline.strip()=="": continue
    glist=gline.split()
    if len(glist) < 9: continue
    if glist[2] != 'gene': continue
    name=list()
    name=glist[8].split(';')
    id=name[0].replace('ID=','')

    try:
        gene_info[glist[0]][id]={'start':int(glist[3]),'end':int(glist[4]),'ID':id}
    except:
        try:
            gene_info[glist[0]][id]={'start':int(glist[3]),'end':int(glist[4]),'ID' : id}

        except:
            gene_info[glist[0]]={id:{'start': int(glist[3]), 'end': int(glist[4]), 'ID' : id }}



#pprint.pprint(gene_info)
#print("Num gff read:",len(gene_info))
#quit()

#### write junction count in a seperate file to be used for DE analysis.
try:
    ofile = ofile
except:
    ofile = "Merged.SJ.table"


##### collect and analyze start and end data of juncions to identify overlapping/splice variants.
sposarray={}
eposarray={}

for chr in gstrt.keys():
    sposarray={chr:[]}
    eposarray = {chr: []}
    spos=0
    epos=0
    for strt in sorted(gstrt[chr].keys()):
        sposarray[chr].extend([strt])


        for end in sorted(gstrt[chr][strt].keys()):
            eposarray[chr].extend([end])
            print(spos, ":", strt, " ", epos, ":", end)
            ### check back 4 steps for overlapping splice variants.
            for i in range(1,5):
                if spos < i or epos < i: continue

                if overlaps([strt,end],[sposarray[chr][spos-i],eposarray[chr][epos-i]]):
                    print(spos,":",strt," ",epos,":",end,"Is splice alter of ",spos-i,":",sposarray[chr][spos-i],epos-i,":",eposarray[chr][epos-i],"*"*i)

        epos += 1
        spos += 1







o = open(ofile, "w")
print("Saving results in :",ofile)

o.write("\t".join(["Chr", "Start", "End", "Start_Gene", "End_Gene",""]))
o.write("\t".join(sorted(sj.keys())))
o.write("\n")
for chr in gstrt.keys():
    for strt in gstrt[chr].keys():
        for end in gstrt[chr][strt].keys():
            sgname = "NoGene"
            egname = "NoGene"
            try:
                for gene in gene_info[chr].keys():
                    if min(int(gene_info[chr][gene]['start']),int(gene_info[chr][gene]['end'])) <= int(strt) <= max(int(gene_info[chr][gene]['start']),int(gene_info[chr][gene]['end'])):
                        sgname=gene
                        #print("%i is in %s start %s - end %s" % (strt,gene,gene_info[chr][gene]['start'],gene_info[chr][gene]['end']) )
                    if min(int(gene_info[chr][gene]['start']),int(gene_info[chr][gene]['end'])) <= int(end) <= max(int(gene_info[chr][gene]['start']),int(gene_info[chr][gene]['end'])):
                        egname=gene
                        #print("%i is in %s start %s - end %s" % (end, gene, gene_info[chr][gene]['start'], gene_info[chr][gene]['end']))
                    if sgname != "NoGene" and egname != "NoGene":
                        break
            except:

                print("NO Gene entry found for %s at loc %i and %i in gff file:" % (chr,int(strt),int(end)))



            o.write("\t".join([chr, str(strt), str(end),sgname,egname]))
            for file in sorted(sj.keys()):
                try:
                    uread = sj[file][chr][strt][end]['uread']
                    #print ("value was found for uread",uread)
                except:
                    uread = 0

                o.write("\t")
                o.write(str(uread))

            o.write("\n")

o.close()
