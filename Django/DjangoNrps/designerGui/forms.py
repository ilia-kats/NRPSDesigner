from designerGui.models import NRP
from django.forms.models import inlineformset_factory, formset_factory
from django.forms import ModelForm, CharField, Textarea
from django.forms.widgets import TextInput


class NRPForm(ModelForm):
	description = CharField(widget=Textarea, required=False)
	class Meta:
		model = NRP
		exclude = ['monomers','owner','designed','construct','designerDomains']