import cookielib
import mechanize
from shutil import copyfile
import sys
import time
import os.path
import os
cwd = os.getcwd()

seq=sys.argv[1]

outfile=cwd+"/"+seq+".mercator_res.gz"
#outfile="/home/ratnesh.singh/Turf/Turf_RNASeq_alternatesplicing/splicegrapher_run/splicegraph_seq_mercator/"+seq+".mercator_res.gz"

## exit if mercator result file exists
if os.path.isfile(outfile):
	sys.exit("Result exists\nMercator result: " + outfile  + " file exists. Please delete or rename before running mercator on file:" + seq + "\n")
elif os.path.isfile(".".join(outfile.split(".")[:-1])):
	sys.exit("Unzipped result file exits\nMercator result: " + ".".join(outfile.split(".")[:-1])  + " file exists. Please delete or rename before running mercator on file:" + seq + "\n")



#sys.exit("\nExiting by choice\n")

##create instance of browser and set cookie handling to cookielib
br=mechanize.Browser()
cj=cookielib.LWPCookieJar()
br.set_cookiejar(cj)

## set browser parameters
br.set_handle_equiv(True)
br.set_handle_gzip(True)
br.set_handle_redirect(True)
br.set_handle_referer(True)
br.set_handle_robots(True)
br.set_handle_robots(False)
br.set_handle_refresh(mechanize._http.HTTPRefreshProcessor(), max_time=1)
br.addheaders=[('User-agent','Mozille/5.0 (X11;Linux i686;en-US;rv:1.9.0.1) Gecko/2008071615 Fedora/3.0.1-1.fc9 Firefox/3.0.1')]


## visit mercator website and select form to fill information
r=br.open('http://www.plabipd.de/portal/web/guest/mercator-sequence-annotation')
br.select_form(nr=0)
## fill required information to form
br.form['Name']='Ratnesh'
br.form['Email']='Ratnesh.Singh@ag.tamu.edu'


## select checkboxes to true
br.find_control('CHLAMY').items[0].selected=True
br.find_control('ORYZA').items[0].selected=True
br.find_control('IPR').items[0].selected=True
br.find_control('IS_DNA').items[0].selected=True
## add file to forms
br.form.add_file(open(seq),'text/plain', seq)

## submit the form 
res=br.submit()

## analyze rsponse
res=br.response()
## to view the response page: print res.read()
## 
while  "Reload to update status" in br.response().read():
	print "Run is not done. Waiting for the response\n"
	time.sleep(10)
	br.select_form(nr=0)
	file=br.submit()
	file=br.response()
	
	
br.select_form(nr=0)
file=br.submit()
file=br.response()	
f=br.retrieve(file.geturl())[0]
copyfile(f,"/home/ratnesh.singh/Turf/Turf_RNASeq_alternatesplicing/splicegrapher_run/splicegraph_seq_mercator/"+seq+".mercator_res.gz")
