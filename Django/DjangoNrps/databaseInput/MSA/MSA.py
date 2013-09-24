import subprocess
import sys
import json

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio import AlignIO




def msa_run(seqRecordList, algorithm = "clustalo"):
	if algorithm == "clustalo":
		# meh write it on our own rather than use BioPython Clustalo wrapper
		# because biopython one does not accept stdin input....
		child = subprocess.Popen("clustalo -i -",
			stdin=subprocess.PIPE,
			stdout=subprocess.PIPE,
			stderr=subprocess.PIPE,
			shell=(sys.platform!="win32"))

		SeqIO.write(seqRecordList, child.stdin, "fasta")
		child.stdin.close()
		clustaloMSA = AlignIO.read(child.stdout, "fasta")
		seqs = msa_to_json(clustaloMSA)
		return seqs

def msa_to_json(alignIO_object):
		msaSeqs = [str(sequence.seq) for sequence in alignIO_object]
		msaNames = [sequence.name for sequence in alignIO_object]
		seqs = json.dumps({'msaSeqs': msaSeqs, 'msaNames': msaNames})
		return seqs