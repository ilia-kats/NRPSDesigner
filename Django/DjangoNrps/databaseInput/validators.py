from django.core.exceptions import ValidationError

def validate_coding_seq(sequence):
	if not len(sequence) % 3:
		msg = "Cds must have length which is a multiple of 3!"
		raise ValidationError(msg)
	if not sequence[0:3].upper() == "ATG":
		msg = "Cds should be starting with ATG!"
		raise ValidationError(msg)

