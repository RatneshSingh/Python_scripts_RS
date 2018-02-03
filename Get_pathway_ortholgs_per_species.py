from bs4 import BeautifulSoup
import mechanize
import re
import urllib2
import pprint
import sys
import os
def get_fasta(link):
    print "Getting fasta from :\t" + link
    dbr=mechanize.Browser()
    dbr.open(link)
    dsoup=BeautifulSoup(dbr.response().read(),'lxml')
    #print dsoup.prettify
    dseq=dsoup.find('pre').contents
    return(dseq[3])


def extract_fasta(burl,href):
    page=urllib2.urlopen(href).read()
    #print "Page:"+page
    nsoup=BeautifulSoup(page,'html.parser')
    for link in  nsoup.find_all('a'):
        glink=link.get('href')
        #print "link.glink:"
        #print glink
        if glink != None and '/dbget-bin/www_bget?-f+-n+a+' in glink:
            print "glink:"+glink
            fasta=get_fasta(str(burl+glink))
            return(fasta)

def get_gene_name(href):
    page=urllib2.urlopen(href).read()
    #print "Page:"+page
    nsoup=BeautifulSoup(page,'html.parser')
    atxt=nsoup.text
    #print atxt
    gmatch=re.search("Gene\s+name[\s\n]+([^\n]+)",atxt)
    gene_name=None
    if gmatch:
        #pprint.pprint(gmatch)
        gene_name=gmatch.group(1)
    else:
        print "***No gene name found"
    return(gene_name)




#htm=open("/home/ratnesh.singh/Turf/Turf_Diamond_Photosynthesis_Pathway/using_kegg/Pathway_00710_plant_species_orthologs.html")
#htm=open("/home/ratnesh.singh/Turf/Turf_Diamond_Photosynthesis_Pathway/using_kegg/Pathway_00195_plant_species_orthologs.html")
#htm=open("/home/ratnesh.singh/Turf/Turf_Diamond_Photosynthesis_Pathway/using_kegg/Pathway_00196_plant_species_orthologs.html")
html=sys.argv[1]
slect=sys.argv[2]
ofolder=sys.argv[3]

slct={}
if slect:
    with open(slect) as sf:
        for st in sf:
            st=st.strip()
            if st:
                slct[str(st).capitalize()]=str(st).capitalize()
    print "Will retrieve sequences for selected species only: " + str(", ".join(slct.keys()))



if not html:
    sys.exit()


htm=open(html)
try:ofolder=sys.argv[3]
except:ofolder=os.getcwd()

burl="http://www.kegg.jp"

soup=BeautifulSoup(htm,'html.parser')

txt=soup.find('thead')

heads=[]

for i in txt.text.split():
    el=i.split('(')
    try:i=el[1]
    except:i=el[0]
    heads.append(str(i.split(')')[0]))


heads=tuple(heads)


trn=0
tdn=0
an=0
for i in soup.find_all('tr'):
    print "\n\nfor tr \n"
    tdn=-1

    for j in i.find_all('td'):
        tdn+=1
        if tdn ==0:grp=j.text
        if tdn ==1:gen=j.text
        if tdn ==2:org=j.text
        if tdn < 3: continue
        if str(org).capitalize() not in slct:
            print "Skipping " + str(org).capitalize() + " as it is not in the list of selected organism"
            continue


        sname=str(str(org).capitalize() + "_" + heads[tdn])  +  " " + str(gen) + " "
        fname=str(heads[tdn])+".prot.fasta"
        try:of=open(ofolder+fname,"a")
        except:of=open(ofolder+fname,"w")
        print "Processing genes for species: " + org
        for k in j.find_all("a"):
            href=k.get('href')
            if 'dbget-bin/www_bget?ko:' in href:continue
            if href != None:
                fasta=extract_fasta(burl,href)
                gname=get_gene_name(href)
                if gname:fasta=fasta.replace(">",">"+str(org).capitalize()+"_"+str(gname)+" ")
                else:fasta=fasta.replace(">",">"+str(sname))
                of.write(fasta)

            else:
                print "href is None:" + href

        of.close()
