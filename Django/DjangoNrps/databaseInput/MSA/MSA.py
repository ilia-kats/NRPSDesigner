import subprocess
import sys
import json

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO

from databaseInput.models import Domain


def clustaloMSA(domType):
	tDomains = Domain.objects.filter(domainType__name = domType)
	okTDomain = [dom for dom in tDomains if len(dom.get_sequence())>0]
	mySeqs = [SeqRecord(Seq(x.get_sequence())) for x in okTDomain]

	# meh write it on our own rather than use BioPython Clustalo wrapper
	# because biopython one does not accept stdin input....
	child = subprocess.Popen("clustalo -i -",
		stdin=subprocess.PIPE,
		stdout=subprocess.PIPE,
		stderr=subprocess.PIPE,
		shell=(sys.platform!="win32"))

	SeqIO.write(mySeqs, child.stdin, "fasta")
	child.stdin.close()
	clustaloMSA = AlignIO.read(child.stdout, "fasta")

	clustaloSeqs = [str(sequence.seq) for sequence in clustaloMSA]
	clustaloNames = [sequence.name for sequence in clustaloMSA]
	seqs = json.dumps({'msaSeqs': clustaloSeqs, 'msaNames': clustaloNames})
	return seqs