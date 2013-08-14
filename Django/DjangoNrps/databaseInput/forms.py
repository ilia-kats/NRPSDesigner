from databaseInput.models import Origin, Cds, Domain
from django.forms.models import inlineformset_factory
from django.forms import ModelForm, HiddenInput


CdsFormSet = inlineformset_factory(Origin, Cds, extra=1)
DomainFormSet = inlineformset_factory(Cds, Domain, extra=5)

class CdsForm(ModelForm):
	class Meta:
		model = Cds


class OriginForm(ModelForm):
	class Meta:
		model = Origin

