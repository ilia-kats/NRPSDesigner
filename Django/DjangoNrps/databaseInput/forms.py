from databaseInput.models import (Origin, Cds, Domain, 
	Substrate, Modification, Product, Linkout, 
	get_pubmed_type, DomainTuple, Experiment)
from django.forms.models import inlineformset_factory, formset_factory
from django.forms import ModelForm, HiddenInput, ModelChoiceField, Form, CharField, Textarea, BooleanField
from django.core.validators import MinLengthValidator

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
	pfamStart = CharField( widget=TextInput(attrs={'class':'disabled', 'readonly':'readonly'}))
	pfamStop = CharField( widget=TextInput(attrs={'class':'disabled', 'readonly':'readonly'}))
	description = CharField(required=False, widget = Textarea(attrs={'rows':5}) , label = "Description:")


	class Meta:
		model = Domain
		fields = ['module', 'domainType', 'substrateSpecificity', 'description',
			'pfamStart', 'pfamStop', 'definedStart','definedStop', 'chirality']

class ProductForm(ModelForm):
	class Meta:
		model = Product
		fields = ['name','description']


# form which enables saving one substrate, or a substrate together with its enantiomer
class DoubleSubstrateForm(ModelForm):
	hasEnantiomer = BooleanField(required=False)
	enantiomerStructure = CharField(required=False , widget = Textarea, label = "Enantiomer Structure:")
	class Meta:
		model = Substrate
		fields = ['name','chirality','structure']

	def __init__(self, user, *args, **kwargs): 
		self.user = user 
		super(DoubleSubstrateForm, self).__init__(*args, **kwargs) 

	def save(self, force_insert=False, force_update=False, commit=True):
		m1 = super(DoubleSubstrateForm, self).save(commit=False)
		m1.user = self.user

		contents = self.cleaned_data


		if contents["hasEnantiomer"]:
			name = m1.name
			m1.name = str(m1.chirality) + "-" + name

			#now instantiate its enantiomer!
			m2args = {'name': m1.reverse_chirality() + "-" + name,
				'chirality': m1.reverse_chirality(),
				'enantiomer':m1,
				'structure': contents["enantiomerStructure"],
				'user': m1.user}

			m2 = Substrate(**m2args)

			m1.enantiomer = m2

			if commit:
				m1.save()
				m2.save()

				m1.enantiomer = m2
				m2.enantiomer = m1

				m1.save()
				m2.save()
		else:
			if commit:
				m1.save()

		return m1

# class to save domain tuple together with its experimental data..
class ExpDomainTupleForm(ModelForm):
	description = CharField(required=True , widget = Textarea, label = "Description:")
	class Meta:
		model = DomainTuple
		exclude = ['experiment']

	def __init__(self, user, *args, **kwargs): 
		self.user = user 
		super(ExpDomainTupleForm, self).__init__(*args, **kwargs) 

	def save(self, force_insert=False, force_update=False, commit=True):
		m1 = super(ExpDomainTupleForm, self).save(commit=False)
		m2 = Experiment(description = self.cleaned_data["description"], user=self.user)
		m1.experiment = m2;
		if commit:
			m2.save()
			# reset experiment of m1 to avoid integrity error
			m1.experiment = m2;
			m1.save()
		return m1


class LinkoutForm(ModelForm):
	identifier = CharField(label="Identifier", validators=[MinLengthValidator(3)])
	
	def __init__(self, *args, **kwargs):
		super(LinkoutForm, self).__init__(*args, **kwargs)
		self.fields['linkoutType'].empty_label = None
		# following line needed to refresh widget copy of choice list
		self.fields['linkoutType'].widget.choices = self.fields['linkoutType'].choices

	def is_valid(self):
		# run the parent validation first
		valid = super(LinkoutForm, self).is_valid()

		if not valid:
			return valid

		tmp = self.save(commit=False)
		try: 
			tmp.linkoutType
		except:
			return False
		return True
	class Meta:
		model = Linkout
		fields = ('linkoutType', 'identifier')


LinkoutFormSet = formset_factory(LinkoutForm)