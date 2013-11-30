from databaseInput.models import Origin, Cds, Domain, Substrate, Modification, Product
from django.forms.models import inlineformset_factory, formset_factory
from django.forms import ModelForm, HiddenInput, ModelChoiceField, Form, CharField, Textarea

from django.forms.widgets import TextInput

CdsFormSet = inlineformset_factory(Origin, Cds, extra=1)


class SubstrateListForm(Form):
	substrate = ModelChoiceField(queryset = Substrate.objects.all())
	
class ModificationsListForm(Form):
    modification = ModelChoiceField(queryset = Modification.objects.all())

SubstrateFormSet = formset_factory(SubstrateListForm, extra=1)
ModificationsFormSet = formset_factory(ModificationsListForm, extra=1)

class CdsForm(ModelForm):
	geneName = CharField(label="Gene name:")
	dnaSequence = CharField(widget = Textarea, label = "DNA sequence:")
	description = CharField(widget = Textarea , label = "Description:")
	class Meta:
		model = Cds

	# clean method now also strips white space of DNA Sequence
	# and also checks if string is actually a FASTA file
	# then copy FASTA file description into entry's description
	def clean(self):
		cleaned_data = super(CdsForm, self).clean()
		if 'dnaSequence' in cleaned_data:
			dnaSequence = cleaned_data.get("dnaSequence")

			if dnaSequence[0] == ">":
				splitDnaSequence= dnaSequence.splitlines()
				description = cleaned_data.get("description","")
				fastaDescription = splitDnaSequence[0][1:]
				cleaned_data["description"] = description + "\nFASTA Description:" + fastaDescription
				dnaSequence = ''.join(splitDnaSequence[1:])
			cleaned_data["dnaSequence"] = ''.join(dnaSequence.split())
		return cleaned_data

class OriginForm(ModelForm):
	class Meta:
		model = Origin
		fields = ['sourceType','source','species','description']

class DomainForm(ModelForm):
    class Meta:
        model = Domain
        fields = ['module', 'domainType', 'substrateSpecificity', 'description',
        	'pfamStart', 'pfamStop', 'definedStart','definedStop']

class ProductForm(ModelForm):
	class Meta:
		model = Product
		fields = ['name','description']