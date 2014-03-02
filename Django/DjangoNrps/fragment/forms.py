from django import forms
from django.conf import settings

########################################## fields

from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA

class SequenceField(forms.CharField):
	def __init__(self, alphabet=IUPACAmbiguousDNA(), **kw):
		self.alphabet = alphabet
		super(SequenceField, self).__init__(**kw)

	def to_python(self, value):
		if value:
			value = value.replace(' ', '')
			return Seq(value, self.alphabet)
		return None

	def clean(self, value):
		value = value.replace(' ', '')
		#clean the chars
		super(SequenceField, self).clean(value)
		#check id the sequence makes sense
		letters = self.alphabet.letters.lower()
		if not letters:
			raise ValueError("Chosen alphabet '%s' does not define any letters." % str(self.alphabet).strip('()'))
		invalid_letters = ''
		for letter in value.lower():
			if letter not in letters:
				if letter not in invalid_letters:
					invalid_letters = invalid_letters + letter
		if len(invalid_letters) > 0:
			s = ''
			if len(invalid_letters) > 1:
				s='s'
			raise forms.ValidationError(u"Sequence Contains invalid letter%s '%s', %s sequences can contain only '%s'" % (s, invalid_letters.upper(), str(self.alphabet).strip('()'), letters.upper()))

		#sequence ok
		return value


############################################# forms

_PART_TYPE_CHOICES = [('', '')]
_PART_TYPE_CHOICES.extend([(x, x) for x in settings.PREDEFINED_PART_TYPES])

class MetaForm(forms.Form):
	name = forms.CharField(	max_length=32,
							widget=forms.widgets.TextInput(attrs={'cols': 32,})
							)
	desc = forms.CharField(	required = False,
							max_length=4096,
							label="Description",
							widget=forms.widgets.Textarea(
								attrs={'cols': 60, 'rows':4,}))

class SequenceForm(forms.Form):
	seq = SequenceField(	max_length=500000,
							label="Sequence",
							widget=forms.widgets.Textarea(
								attrs={'cols': 60,'rows':11,})
							)

class PartTypeForm(forms.Form):
    part_type = forms.ChoiceField(choices=_PART_TYPE_CHOICES, required=False)

class FragmentForm(MetaForm, SequenceForm):
	pass

