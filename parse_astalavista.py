#/usr/bin/python
from sys import argv
import getopt
import sys
import os.path

def code_to_as(code):
    if code.startswith('0,1-2^'):
        astype='skipped_exon'
    elif code.startswith('0,1^2-'):
        astype='retained_introns'
    elif code.startswith('1^,2^'):
        astype='alt_donors'
    elif code.startswith('1-2^,3-4^'):
        astype='mutually_exclusive_exons'
    elif code == '1-,2-':# or code.startswith('1-,2-'):
        astype='alt_acceptors'
    elif code.startswith('1['):
        astype='alt_TSS'
    else:
        astype='others'

    return(astype)






usage='''
usage:
python3 {} -g ASTALVISTA_gtf_files_seperated_with_comma   -o outfile_to_save_results.

'''.format(argv[0])

try:
    optlist,arg=getopt.getopt(argv[1:],'g:o:h')
except:
    print("Incorrect commandline options "+ argv[1:])
    print(usage)

## initialize variables
out=None
files=None

for o, a  in optlist:
    if o == '-g':
        files = a
    elif o == '-o':
        out= a
    elif o == '-h':
        print(usage)
        sys.exit()
    else:
        assert False, "Unhandle option"
        print(usage)

## check opetions
if not files:
    print(usage)
    exit()

#of=open(gtf,r)
asta=dict()
askeys=dict()
totaskeys=dict()
astypes=dict()
ast_list=dict()
filelist=[]



## read gtf file and collect relevant info in dict.
for gtf in files.split(","):
    if gtf.strip() == "":continue
    filelist.append(gtf.strip())
    print("\nProcessing file: {}".format(gtf))
    with open(gtf,'r') as gf:
        asta[gtf]={}
        astypes[gtf]={}
        totaskeys[gtf]=0
        for line in gf:
            line=line.strip()
            col=line.split("\t")
            prpt=col[8].split(";")
            for s in prpt:
                if not s.strip() : continue
                k,v,ext = s.split('"')
                if  "structure" in k:
                    askeys[v.strip()]=1
                    totaskeys[gtf]+=1
                    try: asta[gtf][v.strip()]+=1
                    except:asta[gtf][v.strip()]=1
                    asty=code_to_as(v.strip())
                    try: astypes[gtf][asty]+=1
                    except:astypes[gtf][asty]=1
                    ast_list[asty]=asty

                    #print("Processing propt elements: " + " :" + str(asta[gtf][v.strip()]) + ":" + v.strip())

if out:
    outn=out.rsplit(".")
    if len(outn)>1:
        out1=".".join(out.rsplit(".")[0:-1])+".abs."+out.rsplit(".")[-1]
        out2=".".join(out.rsplit(".")[0:-1])+".percent."+out.rsplit(".")[-1]
        out3=".".join(out.rsplit(".")[0:-1])+".AStypeSummary.abs."+out.rsplit(".")[-1]
        out4=".".join(out.rsplit(".")[0:-1])+".AStypeSummary.percent."+out.rsplit(".")[-1]
        i=1
        while os.path.isfile(out1):
            out1=".".join(out.rsplit(".")[0:-1])+".abs." + str(i) + "." + out.rsplit(".")[-1]
            out2=".".join(out.rsplit(".")[0:-1])+".percent." + str(i) + "." + out.rsplit(".")[-1]
            out3=".".join(out.rsplit(".")[0:-1])+".AStypeSummary.abs." + str(i) + "." + out.rsplit(".")[-1]
            out4=".".join(out.rsplit(".")[0:-1])+".AStypeSummary.percent." + str(i) + "." + out.rsplit(".")[-1]
            i+=1
    else:
        out1=out+".abs.table"
        out2=out+".percent.table"
        out3=out+".AStypeSummary.abs.table"
        out4=out+".AStypeSummary.percent.table"
        i=1
        while os.path.isfile(out1):
            out1=out+".abs." + str(i) + "." + out.rsplit(".")[-1]
            out2=out+".percent." + str(i) + "." + out.rsplit(".")[-1]
            out3=out+".AStypeSummary." + str(i) + ".abs." + out.rsplit(".")[-1]
            out4=out+".AStypeSummary." + str(i) + ".percent." + out.rsplit(".")[-1]
            i+=1

else:
    out1="Astalavista_event_summary.abs.table"
    out2="Astalavista_event_summary.percent.table"
    out3="Astalavista_event_summary.AStypeSummary.abs.table"
    out4="Astalavista_event_summary.AStypeSummary.percent.table"
    i=1
    while os.path.isfile(out1):
        out1="Astalavista_event_summary.abs." + str(i) + ".table"
        out2="Astalavista_event_summary.percent." + str(i) + ".table"
        out3="Astalavista_event_summary.AStypeSummary.abs." + str(i) + ".table"
        out4="Astalavista_event_summary.AStypeSummary.percent." + str(i) + ".table"
        i+=1

of1=open(out1,"w")
of2=open(out2,"w")
of3=open(out3,"w")
of4=open(out4,"w")

of1.write("AS_event\t")
of1.write("\t".join(filelist))
of1.write("\nTotal_num_as\t")
of1.write("\t".join([str(totaskeys[x]) for x in filelist]))


of2.write("AS_event\t")
of2.write("\t".join(filelist))
of2.write("\nTotal_num_as\t")
of2.write("\t".join([str(totaskeys[x]) for x in filelist]))
#for file in filelist:
#    of.write("\t{}".format(file))

for events in askeys.keys():
    of1.write("\n{}".format(events.strip()))
    of2.write("\n{}".format(events.strip()))
    for file in filelist:
        try:of1.write("\t{}".format(asta[file][events]))
        except:of1.write("\t{}".format(0))
        try:of2.write("\t{}".format(round(100*asta[file][events]/totaskeys[file],4)))
        except:of2.write("\t{}".format(0))

of1.close()
of2.close()


#print astype Summary
ascode_list=[x for x in ast_list.keys()]
print("File\t"+"\t".join(ascode_list))
of3.writelines("File\t"+"\t".join(ascode_list))
of4.writelines("File\t"+"\t".join(ascode_list))
for file in filelist:
    print("\n{}".format(file),end="")
    of3.writelines("\n{}".format(file))
    of4.writelines("\n{}".format(file))
    for ast in ascode_list:
        print("\t{}({}%)".format(astypes[file][ast],round(100*astypes[file][ast]/totaskeys[file],2)),end="")
        of3.writelines("\t{}".format(astypes[file][ast]))
        of4.writelines("\t{}".format(round(100*astypes[file][ast]/totaskeys[file],2)))

of3.close()
of4.close()
print("\nResults are saved in files:\nAbsolute values:{}\nPercent values:{}\nSummary_abs:{}\nSummary_percent:{}\n\n".format(out1,out2,out3,out4))
