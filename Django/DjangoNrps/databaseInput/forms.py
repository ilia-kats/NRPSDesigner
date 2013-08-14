from databaseInput.models import Origin, Cds, Domain, Substrate
from django.forms.models import inlineformset_factory, formset_factory
from django.forms import ModelForm, HiddenInput, ModelChoiceField, Form

from django.forms.widgets import TextInput

CdsFormSet = inlineformset_factory(Origin, Cds, extra=1)
DomainFormSet = inlineformset_factory(Cds, Domain, extra=5)


class SubstrateListForm(Form):
	substrate = ModelChoiceField(queryset = Substrate.objects.all())

SubstrateFormSet =formset_factory(SubstrateListForm, extra=1)

class CdsForm(ModelForm):
	class Meta:
		model = Cds


class OriginForm(ModelForm):
	class Meta:
		model = Origin
