import requests
import time
from xml.dom.minidom import parseString
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#first read FASTA file and translate sequence to protein!
dnaSeq  = SeqIO.read("bpsa.fasta", "fasta",IUPAC.unambiguous_dna).seq
protSeq = dnaSeq.translate(to_stop=True)

#send PFAM request 1
pfamUrl = "http://pfam.sanger.ac.uk/search/sequence"
pfamParams = {'seq':str(protSeq), 'output':'xml'}
pfamRequest= requests.get(pfamUrl, params=pfamParams)

#extract result link from XML by converting to DOM, extracting result url tag, removing tag elements from string
pfamDom = parseString(pfamRequest.text)
pfamResultUrl = pfamDom.getElementsByTagName('result_url')[0].toxml()
pfamResultUrl = pfamResultUrl.replace("<result_url>","").replace("</result_url>","")

#wait a bit
time.sleep(1)

# keep sending PFAM request 2 until it works
while True:
	pfamResultRequest = requests.get(pfamResultUrl)
	if pfamResultRequest.status_code == 200:
		break
	elif pfamResultRequest.status_code == 202:
		time.sleep(1)
	else:
		break #THIS SHOULD actually throw exception!!!

# read xml DOM, find stuff corresponding to domains
pfamResultDom = parseString(pfamResultRequest.text)
pfamResultMatches = pfamResultDom.getElementsByTagName('match')
pfamDomains = []
for match in pfamResultMatches:
	pfamDomains.append(match.attributes["id"].value)




