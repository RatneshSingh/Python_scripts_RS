import cookielib
import mechanize
from shutil import copyfile
import sys
import time
import os.path
import os
cwd = os.getcwd()

seq=sys.argv[1]
sequence={}
header=""
with open(seq) as F:
	for line in F:
		if line.startswith(">"):header=line.split()[0].replace(">","")
		else:
			try:sequence[header]=sequence[header]+line.strip()
			except:sequence[header]=line.strip()

outfile=""
outprimer=""

#sys.exit("\nExiting by choice\n")

for head in sequence:
	print "running for {} with sequence\n{}".format(head,sequence[head])
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


	## visit primer3 website and select form to fill information
	#url="http://www.bioinformatics.nl/cgi-bin/primer3plus/primer3plus.cgi"
	#url="http://bioinfo.ut.ee/cgi-bin/primer3-0.4.0/primer3_results.cgi"
	url="http://bioinfo.ut.ee/primer3/"
	size="150-250"
	#sequence="ATGGACCGCGGCGAGGTCTCCCGGATGGCGAGCCTCCGTCGGGAGAGCTCGCTGTGGCGGCGCGGCGACGACGTGTTCTGCGGTCGTCGTCGCGGTTCCaggacgaggaggacgacgaggaggCGCTGCGCTGGGCGGCGCTGGAGCGGCTGCCGACCTACGACCGCGTGCGGCGCGGGATCCTGTCCCTGGACGAGGACGGCGAGAAGGTGGAGGTGGACGTCGGCCGCCTCGGAGCGCGCGAGTCGCGCGCGCTCATCGAGCGACTCGTCCGCGCCGCCGACGACGACCACGAGCGGTTCCTGCTCAAGCTCAAGGAGCCATCATGAAGTTCAACTTCCAGAAGAGATGA"

	r=br.open(url)
	br.select_form(id="primer3web")
	## fill required information to form
	br.form['SEQUENCE_ID']=seq+"_"+head
	br.form['SEQUENCE_TEMPLATE']=sequence[head]
	br.form['PRIMER_PRODUCT_SIZE_RANGE']=size


	## submit the form 
	res=br.submit()

	## analyze rsponse
	res=br.response()
	data=res.get_data()

	from BeautifulSoup import BeautifulSoup

	soup = BeautifulSoup(data)
	cols = soup.findAll('pre')

	data=cols[0].renderContents() # print content of first <td> element
	RP=""
	LP=""
	for i in data.split("\n"):
		if "LEFT PRIMER" in i:
			LP=i.split()[-1]
		elif "RIGHT PRIMER" in i:
			RP=i.split()[-1]

		if "ADDITIONAL OLIGOS" in i:
			break


	outfile=cwd+"/"+seq+head+".primer3_out.txt"
	outprimer=cwd+"/"+seq+head+".primer3.txt"
	with open(outprimer,"w") as OP, open(outfile,"w") as OF :
		print OP.write("{}\t{}_F\t{}\n{}\t{}_R\t{}\n".format(seq,head,LP,seq,head,RP))
		print OF.write(data)


