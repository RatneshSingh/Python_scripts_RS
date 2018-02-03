#!/usr/bin/python
from sys import argv, exit
from BCBio import GFF
#from collections import defaultdict
# from BCBio.GFF import GFFExaminer
import pprint
import getopt

### functions
def usage():
    usages='''
python {} options....
-h  print help and exit.
-f  table files seperated with comma.
-o  outfile to save merged result
-s  keep the same order of files as provided.
-k  key column to use for comparison [1]
-c  column to merge [All]
-m  missing value replaced with [0]
-d  Split columns with delimiter

'''.format(argv[0])
    print(usages)
    exit()






try:
    opts, args = getopt.getopt(argv[1:], 'hg:f:o:m:k:d:c:s')
except:
    print(argv[0], " -f infile -o  output_file -k col_to_use_as_key_to_join[1]  -c col_to_join[All]  -m show_missing_values_as[0]  -d delim_to_parse[\s]")
    print("One of the argument in:", argv, " is invalid:")
    exit()
files = str()
mmis=0
mkey=1
mcol="All"
keep_order=False

for opt, arg in opts:
    if opt == "-h":
        usage()
        exit(" ")
    elif opt in ("-f"):
        files = arg.split(',')
        #print("Found files:", files)
    elif opt in ("-o"):
        ofile = arg
        print("Results will be saved in", ofile)
    elif opt in ("-k"):
        mkey = int(arg)
        #print("Files will be merged based on column:", mkey)
    elif opt in ("-c"):
        mcol = int(arg)
        print("The column being merged:", mcol)
    elif opt in ("-m"):
        mmis = arg
        #print("Replace missing values with:", mmis)
    elif opt in ("-d"):
        mdim = arg
        print("Split columns with delimiter:", mdim)
    elif opt in ("-s"):
        keep_order=True
        print("Output will have files contents in same order as input order.")
    else:
        print("Use legal flags\n", argv)
        print("error:", args)
        exit("Unable to parse argument")

if len(files) < 1:
    print("No input file found to process")
    quit()


if mkey > 0: mkey=mkey - 1
else: print "-k values can only be positive integer. Defaulting to k=1"



print "Col number "+str(mcol)+" from each files will be merged using keys in col :"+str(mkey+1)+"\nMissing values will be replaced by:"+str(mmis)





sj={}
flist=[]
klist={}
smis={}
file_in_order=[]
## read file to collect sj boundry and read count
for file in files:
    if file.strip()=="":continue
    sj[file]={}
    file_in_order.append(file);
    smis[file]={}
    flist.append(file)
    print("Reading ",file)
    with open(file,"r") as f:
        for line in f:
            elem=[]
            elem=line.strip().split()
            if isinstance(mcol,int) and mcol-1 >= len(elem):
                print("The number of columns ({}) in file is less than requested number {}".format(len(elem),mcol))
                quit()
            if mkey >= len(elem):
                print("The number of columns ({}) in file is less than requested key {}".format(len(elem),mkey+1))
                quit()

            klist[elem[mkey]]=1
            if isinstance(mcol,int):
                sj[file][elem[mkey]]=elem[mcol-1]
                smis[file]['mmis']=str(mmis)
            else:
                sj[file][elem[mkey]]="\t".join([x for x in elem if elem.index(x) != mkey ])
                smis[file]['mmis']=(str(mmis)+"\t") * ( len(elem) -1 )


#pprint.pprint(sj)
#quit()
print("Number of files read:",str(len(sj)))
#pprint.pprint(sj)
#quit()


#### write merged columns in a seperate file to be used for DE analysis.
try:
    ofile = ofile
except:
    ofile = "Merged."+str(mcol)+"columns_by_keyCol_"+str(mkey)+".table"

o = open(ofile, "w")
print("Saving results in :",ofile)

if not keep_order: file_in_order=sorted(sj.keys())

o.write("\t")
o.write("\t".join(file_in_order))
o.write("\n")

for tkey in sorted(klist.keys()):
    o.write(tkey)
    o.write("\t")
    for tfile in file_in_order:
        try:
            o.write(sj[tfile][tkey])
        except:
            o.write(smis[tfile]['mmis'])
        o.write("\t")
    o.write("\n")

o.close()
