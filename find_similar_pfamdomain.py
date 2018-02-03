import argparse
import sys
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--query', help='domtblout file for query sequence obtained from hmmscan run.')
group.add_argument('--qname', help='query name to be used as ref from domtblout file. only in case seperate tblout for query file is not provided and query is included in the homolog tblout')
parser.add_argument('--homologs','--hmlg', required=True, help='domtblout file for homologous sequence obtained from hmmscan run.')
parser.add_argument('--out', help='outfile to save list of homologous sequences qualify as homolog')
parser.add_argument('--ival', help='independent e-value cutoff to assign a domian to a sequence',type=float, default=1e-5)
parser.add_argument('--perc_dom','-pd',help='Percent of query domains to be shared by a sequence to qualify as a homolog', type=int,default=80)
args = parser.parse_args()

query={}
hmlg={}
out=""
if args.out is not None:
    out=args.out
else:
    out="hmmscan_homologs_of_"+args.homologs+".list"


print("{}")
### read domtblout for query sequences.
if (args.query is not None):
    with open(args.query) as Q:
        for line in Q:
            if line.startswith("#"): continue
            elem=line.split()
            if float(elem[12]) <= float(args.ival):
                try:query[elem[3]].update({elem[1]:elem[12]})
                except:query[elem[3]]={elem[1]:elem[12]}

### read domtblout for homologous sequences.
with open(args.homologs) as Q:
    for line in Q:
        if line.startswith("#"): continue

        elem=line.split()
        if (args.qname is not None and elem[3].strip() == args.qname.strip()):
            try:query[elem[3]].update({elem[1]:elem[12]})
            except:query[elem[3]]={elem[1]:elem[12]}

        elif float(elem[12]) <= args.ival:
            try:hmlg[elem[3]].update({elem[1]:elem[12]})
            except:hmlg[elem[3]]={elem[1]:elem[12]}


if len(query) < 1:
    print("No query sequence found. Exiting.")
    sys.exit()

## Process domians.

selected={}
for qry in query:
    print("\n\n###################################")
    print(qry,query[qry])
    for seq in hmlg:
        dom_tot=0
        dom_in=0
        for doms in query[qry]:
            dom_tot+=1
            if doms in hmlg[seq]:
                dom_in+=1

        if dom_in*100/dom_tot >= args.perc_dom:
            print("\n{} \tshared \t{} \tdomains with \t{} at ivalue <= {}".format(seq,dom_in,qry,args.ival))
            selected[seq]=1

with open(out,"w") as OF:
    for sname in selected:
        OF.write("{}\n".format(sname))
