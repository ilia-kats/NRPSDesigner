from django.core.exceptions import ValidationError
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Data.CodonTable import TranslationError

def validateCodingSeq(sequence):
	try:
		Seq(sequence, IUPAC.unambiguous_dna).translate(table="Bacterial", cds=True)
	except TranslationError as e:
		raise ValidationError(e)
