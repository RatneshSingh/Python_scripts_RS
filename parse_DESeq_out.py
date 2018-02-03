#!/usr/bin/python
from sys import argv, exit
from BCBio import GFF
from Bio import SeqIO
#from collections import defaultdict
# from BCBio.GFF import GFFExaminer
import pprint
import getopt
import copy
import sys
### functions

def overlaps(frst,scnd):
    frst = sorted(frst)
    scnd = sorted(scnd)
    return min(frst[1],scnd[1]) - max(frst[0],scnd[0]) > 0

def usage():
    help="""\n ...options....
    
    -d  deseq_result_file containing results from DESeq2.
        Result file contains columns: BaseMean, log2foldchange, lfcSE, stat, pvalue, padj, Sample_Names.....
    -j  Junction_file. Merged SJ junction file containing:
        Chr,start,end, Start_Gene,End_Gene and NumReads for each sample with sample name as column header.
        This file could be created by using merge_START_SJout.py script
    -o  output_file to save resulting table result.[auto generated]
    -t  padj|pval. Use this for filtering [ padj]
    -p  pvalue_limit [0.05]. Dont report junctions below this significance level
    -l  log2_foldchange_difference between junctions [2]
    -m  max_dist_between_junction. junctions farther than this will not be reported [10000].
        This is to avoid in-correct mapping causing invalid junctions to be reported.
    -f  fasta_file_to_extract_seqs. To extract sequence containing the junctions for further analysis.
    -b  depth of loopback to search for overlapping junctions [5]
        every junction is compared with this many previous junctions to find overlap.
    -v  print progress report.
    -h  Print help and exit.
    
    """
    sys.stderr.write("\n"+argv[0] + " -d deseq_result_file|-j Junction_file  -o  output_file -t padj|pval -p  pvalue_limit [0.05] -l log2_foldchange_difference[2] -m max_dist_between_junction[10000] -f fasta_file_to_extract_seqs \n"+help)
    exit()

def get_seq_region():
    print("Read fasta sequences into dictionary to ")
    if seq in dGeno:
        return(seq,dGeno[seq])
    else:
        print("This part is still under construction\n")
        

def parse_DESeq(defile):
    ginfo={}
    header=None
    with open(defile,"r") as df:
        line=str()
        for line in df:
            line=line.strip()
            lel=line.split()
            names=lel[0].strip('"').split(":")
            #pprint.pprint(names)
            cont=names[0]
            if cont == "baseMean":
                header="\t".join(lel)
                header="Gene_Junc\t"+header
                continue
            strt=int(names[1])
            end=int(names[2])
            gst=names[3]
            gnd=names[4]
            log2fc=float(lel[2])

            try:pval=float(lel[5])
            except:pval=lel[5]

            try: padj = float(lel[6])
            except: padj = lel[6]

            try:
                ginfo[cont][strt].update({end:{'gst':gst,'gnd':gnd,'data':line,'log2fc':log2fc,'pval':pval,'padj':padj}})
            except:
                try: ginfo[cont].update({strt:{end:{'gst':gst,'gnd':gnd,'data':line,'log2fc':log2fc,'pval':pval,'padj':padj}}})
                except:  ginfo.update({cont:{strt:{end:{'gst':gst,'gnd':gnd,'data':line,'log2fc':log2fc,'pval':pval,'padj':padj}}}})
    return(ginfo,header)
    df.close()

## Parse merged junction file containing Chr,start,end and NumReads for each sample with sample name as column header.
def parse_junct(jefile):
    jinfo={}
    header=[]
    with open(jefile,"r") as jf:
        for line in jf:
            line=line.strip()
            names=line.split()
            cont=names[0]
            if cont == "Chr" and names[1] == "Start" and names[2] == "End":
                header=[x.strip() for x in names[5:]]
                continue
            strt=int(names[1])
            end=int(names[2])
            gst=names[3]
            gnd=names[4]
            cnt=[int(x) for x in names[5:]]
            for i in range(0,len(cnt)):
                try:
                    jinfo[header[i]][cont][strt].update({end:{'gst':gst,'gnd':gnd,'data':cnt[i]}})
                    #print("\t\tAdding Start {}:key {} in hash table".format(i,strt))
                except:
                    try:
                        jinfo[header[i]][cont].update({strt:{end:{'gst':gst,'gnd':gnd,'data':cnt[i]}}})
                        #print("\tStarting for sample {} new contig {}: as key {} in hash table".format(header[i],cont,i))
                    except:
                        try:
                            jinfo[header[i]].update({cont:{strt:{end:{'gst':gst,'gnd':gnd,'data':cnt[i]}}}})
                        except:
                            jinfo.update({header[i]:{cont:{strt:{end:{'gst':gst,'gnd':gnd,'data':cnt[i]}}}}})
                        #print("Starting new {}:key {} in hash table".format(i,header[i]))
    jf.close()
    print("Read {} contigs in file {}".format(len(jinfo),jefile))
    return(jinfo,header)
    

def get_overlapping_junc(sdb,stpos,tst1,ten1,sbreadth,smdist):
    stpos=int(stpos)
    tst1=int(tst1)
    ten1=int(ten1)
    sbreadth=int(sbreadth)
    ov_list=[]
    index_list=[]
    ## find overlap from current strt,end pair to previous sbreadth pairs
    for i in range(0,sbreadth+1):
        if stpos < i:continue
        for j in range(0,len(sdb[stpos -i][1])):
            tst2=sdb[stpos -i][0]
            ten2=sdb[stpos -i][1][j]
            if tst1 == tst2 and ten1 == ten2: continue
            if abs(tst2 -ten2) >= smdist: continue
            if overlaps([strt,end],[tst2,ten2]):
                ov_list.append([tst2,ten2])
                index_list.append([stpos -i,j])
                if verbose:print("Found overlap between {} - {}  and  {} - {}".format(tst2,ten2, tst1 ,ten1))
        
    return(ov_list,index_list)            



####################################################################################

try:
    opts, args = getopt.getopt(argv[1:], 'ho:d:p:b:l:t:m:f:j:v')
except:
    sys.stderr.write("\nOne of the argument in:"+ argv +" is invalid:")
    usage()
### Defaults
files = str()
diff_lfc = 2
lim_pval = 0.05
cnt_lim = 1
pstr = 'padj'
mdist = 10000
breadth=5
ofile=None
dfile=None
jfile=None
fasta_file=None
geno=None
verbose=None
### parse commandline options
for opt, arg in opts:
    if opt == "-h":
        usage()
    elif opt in ("-o"):
        ofile = arg
        sys.stderr.write("\n***Results will be saved in " +  ofile)
    elif opt in ("-d"):
        dfile = arg
        sys.stderr.write("\n***Will read DESeq results from " + dfile)
    elif opt in ("-j"):
        jfile = arg
        sys.stderr.write("\n***Will read Junction information from " + jfile)
    elif opt in ("-l"):
        diff_lfc = int(arg)
        sys.stderr.write("\n***log2fold change cutoff >= " + str(diff_lfc))
    elif opt in ("-p"):
        lim_pval = float(arg)
    elif opt in ("-t"):
        pstr = arg
    elif opt in ("-v"):
        verbose = 1
    elif opt in ("-f"):
        fasta_file = arg
    elif opt in ("-b"):
        breadth = arg
    elif opt in ("-m"):
        mdist = int(arg)
        sys.stderr.write("\n***Junctions farther than " + str(mdist) +" will not be printed")
    else:
        sys.stderr.write("\nUse legal flags\n" + argv)
        sys.stderr.write("\nerror:" +  args)
        exit("\nUnable to parse argument")


### check if input file is provided
if not dfile and not jfile:
    print("Atleas one input file is required;")
    usage()

if pstr == 'padj' or pstr == 'pval':
    pass
else:
    sys.stderr.write("\n-p can only take values padj or pval. Restoring deafualt value:padj")
    pstr='padj'

sys.stderr.write("\n***will print for " + pstr + "<=" + str(lim_pval) + "max distance between junctions:" + str(mdist))


###
printed={}
sj={}

## read deseq out file

## read merged junction file

### read fasta file
if fasta_file: dGeno=SeqIO.to_dict(fasta_file,"fasta")



### process read data and write results.
## create output file names if not given by user.

if ofile:
    dofile = ofile.strip()
    jofile = ofile.strip()
    if dfile and jfile:
        dofile=str(ofile).rstrip(".table") + ".DE_lfc" + str(diff_lfc) + ".table"
        jofile=str(ofile).rstrip(".table") + ".Junction_summary" + ".table"
else:
    if dfile: dofile=str(dfile).rstrip(".table") + ".DE_lfc" + str(diff_lfc) + ".table"
    if jfile: jofile=str(jfile).rstrip(".table") + ".Junction_summary" + ".table"
    #sys.stderr.write("\n***Results will be saved in " + ofile + "\n")




if dfile:
    of=open(dofile,"w")
    (ginfo,header)=parse_DESeq(dfile)
    print("\nJunc:{:>5}:{:>50}:{:>19}:{:>9}:{:>9}:{:>9}: <====> Junc:{:>5}:{:>9}:{:>9}:{:>10}:".format("Junc1","Gene","Chr","Start","end","Log2FC","Junc2","Start","end","Log2FC"))
    of.write("Junction1\tGene\tChr\tStart\tend\tLog2FoldChange\tJunction2\tStart\tend\tLog2FoldChange\n")

    ef=open(dofile + ".expData.table","w")
    ef.write(header)

    for cont in ginfo.keys():
        scount=-1
        ecount=-1
        scount_db=[]
        ecount_db = []
        for strt in sorted(ginfo[cont].keys()):
            if isinstance(strt,str):continue
            scount_db.append([strt,[]])
            scount += 1
            for end in sorted(ginfo[cont][strt].keys()):
                #print("Scount:{}  ".format(scount))
                scount_db[scount][1].append(end)
                if isinstance(ginfo[cont][strt][end][pstr],str):continue
                if ginfo[cont][strt][end][pstr] > lim_pval: continue
                ov_list,ind_list=get_overlapping_junc(scount_db,scount,strt,end,breadth,mdist)
                if len(ov_list) < 1:continue
                if verbose:pprint.pprint(ov_list)
                for i in range(0,len(ov_list)):
                    ovel=ov_list[i]
                    lfc1=ginfo[cont][strt][end]['log2fc']
                    lfc2=ginfo[cont][ovel[0]][ovel[1]]['log2fc']
                    gst1=ginfo[cont][strt][end]['gst']
                    gnd1=ginfo[cont][strt][end]['gnd']
                    gst2=ginfo[cont][ovel[0]][ovel[1]]['gst']
                    gnd2=ginfo[cont][ovel[0]][ovel[1]]['gnd']
                    dline1=ginfo[cont][strt][end]['data']
                    dline2=ginfo[cont][ovel[0]][ovel[1]]['data']
                    if abs(lfc1 - lfc2) >= diff_lfc:
                        if abs(strt - end) >= mdist: continue
                        if abs(ovel[0] - ovel[1]) >= mdist: continue
                        #print("Entry:",scount,":",cont,scount_db[scount],"lfc:",lfc1," <====>","Entry:",scount-i,":",scount_db[scount-i]," lfc:",lfc2)
                        print("Junc:{:5}:{:>50}:{:>20}{:10}{:10}{:>10.4f} <====> Junc:{:5}:{:10}{:10}{:>10.4f}".format(scount,gst1,cont,strt,end,lfc1,ind_list[i][0],ovel[0],ovel[1],lfc2))
                        of.write("\n{}\t{}\t{}\t{}\t{}\t{:.3f}\t{}\t{}\t{}\t{:.3f}".format(scount,gst1,cont,strt,end,lfc1,ind_list[i][0],ovel[0],ovel[1],lfc2))
                        if dline1 in printed:
                            printed[dline1] += 1
                        else:
                            ef.write("\n"+dline1)
                            printed[dline1] = 1
                        if dline2 in printed:
                            printed[dline1] += 1
                        else:
                            ef.write("\n"+dline2)
                            printed[dline2] = 1
                
    of.close()
    ef.close()

###
## create output file names if not given by user.
if jfile:
    of=open(jofile,"w")
    junc_exp={}
    junc_num_exp={}

    junc_splice={}
    junc_num_splice={}

    junc_gene={}
    junc_num_gene={}

    (jinfo,header)=parse_junct(jfile)
    print("\nTotal number of samples read: {}\t".format(len(jinfo)))
    of.write("SampleName\tTotal_Junctions\tJunctionsWithSplicing\tPercentSpliced\tSplicedGenes")

    #ef=open(ofile + ".expData.table","w")
    #ef.write(header)
    for samp in jinfo.keys():
        if verbose: print("\nProcessing sample:{}\t total contigs:{} contig Names: {}".format(samp,len(jinfo[samp]),"*".join([x for x in jinfo[samp].keys()][:50])))

        junc_exp.update({samp:{}})
        junc_num_exp.update({samp:0})

        junc_splice.update({samp:{}})
        junc_num_splice.update({samp:0})

        junc_gene.update({samp:{}})
        junc_num_gene.update({samp:0})


        for cont in jinfo[samp].keys():
            scount=-1
            scount_db=[]
            ecount_db = []
            
            print("\tProcessing Sequence:{}\t with total junction start sites:{}".format(cont,len(jinfo[samp][cont])))
            for strt in sorted(jinfo[samp][cont].keys()):
                scount_db.append([strt,[]])
                scount += 1
                
                for end in sorted(jinfo[samp][cont][strt].keys()):
                    
                    scount_db[scount][1].append(end)
                    cnt1=jinfo[samp][cont][strt][end]['data']
                    if cnt1 < cnt_lim:
                        if verbose: print("Count1: {} is smaller than limit : {}.. Skipping".format(cnt1,cnt_lim))
                        continue

                    if abs(strt - end) >= mdist:
                        if verbose: print("Distance between junc start {} and end {} is less than {} .. Skipping".format(strt,end,cnt1))
                        continue
                    gst1=jinfo[samp][cont][strt][end]['gst']
                    gnd1=jinfo[samp][cont][strt][end]['gnd']
                    
                    junc_exp[samp][gst1] = 1
                    junc_num_exp[samp] += 1
                    ov_list,ind_list=get_overlapping_junc(scount_db,scount,strt,end,breadth,mdist)
                    if len(ov_list) < 1:continue
                    if verbose: pprint.pprint(ov_list)
                    for ovel in ov_list:
                        cnt2=jinfo[samp][cont][ovel[0]][ovel[1]]['data']
                        if cnt2 < cnt_lim: continue
                        gst2=jinfo[samp][cont][ovel[0]][ovel[1]]['gst']
                        gnd2=jinfo[samp][cont][ovel[0]][ovel[1]]['gnd']
                        
                        junc_splice[samp][gst2]=1
                        junc_num_splice[samp]+=1
                        
                        junc_gene[samp].update({gst1:1})
                        junc_gene[samp].update({gst2:1})
                        
                        if verbose: print("Found overlap between {} - {}  and  {} - {}".format(ovel[0],ovel[1], strt ,end)) 
        print("\n{}\tTotal_Junc:  {}\tSplice_junc:  {}  {:.2f}%".format(samp,junc_num_exp[samp],junc_num_splice[samp],junc_num_splice[samp]*100/junc_num_exp[samp]))
        of.write("\n{}\t{}\t{}\t{:.2f}\t{}".format(samp,junc_num_exp[samp],junc_num_splice[samp],junc_num_splice[samp]*100/junc_num_exp[samp],len(junc_gene[samp])))
        #pprint.pprint(junc_exp[samp])
    #of.close()
    #ef.close()


#pprint.pprint(jinfo)
