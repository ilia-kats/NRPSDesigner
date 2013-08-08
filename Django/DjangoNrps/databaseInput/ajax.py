import simplejson
import requests
import time
import pdb
from xml.dom.minidom import parseString
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from dajaxice.decorators import dajaxice_register

@dajaxice_register
def args_example(request, text):
	#return simplejson.dumps({'message':'Your message is %s!' % text})
	return text

@dajaxice_register
def sauceFunc(request, sequence):
#first read FASTA file and translate sequence to protein!
#dnaSeq  = SeqIO.read("bpsa.fasta", "fasta",IUPAC.unambiguous_dna).seq
	sequence = sequence.replace("\n","").replace("\t","")
	dnaSeq = Seq(sequence, IUPAC.unambiguous_dna)
	protSeq = dnaSeq.translate(to_stop=True)
	
#send PFAM request 1
	pfamUrl = "http://pfam.sanger.ac.uk/search/sequence"
	pfamParams = {'seq':str(protSeq), 'output':'xml'}
	pfamRequest= requests.get(pfamUrl, params=pfamParams)

	#extract result link from XML by converting to DOM, extracting result url tag, removing tag elements from string
	pfamDom = parseString(pfamRequest.text)
	pfamJobId = pfamDom.getElementsByTagName('job')[0].attributes["job_id"].value
	pfamResultUrl = "http://pfam.sanger.ac.uk/search/sequence/resultset/" + pfamJobId
	pfamGraphicUrl = "http://pfam.sanger.ac.uk/search/sequence/graphic/" + pfamJobId

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
	#pfamResultDom = parseString(pfamResultRequest.text)
	#pfamResultMatches = pfamResultDom.getElementsByTagName('match')
	#pfamDomains = []
	#for match in pfamResultMatches:
	#	pfamDomains.append(match.attributes["id"].value)
	pfamResultRequest = requests.get(pfamGraphicUrl)
	#pdb.set_trace()
	return pfamResultRequest.text[1:-1]
	