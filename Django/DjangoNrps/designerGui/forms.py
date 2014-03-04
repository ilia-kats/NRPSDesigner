from designerGui.models import NRP
from django.forms.models import inlineformset_factory, formset_factory
from django.forms import ModelForm, CharField, Textarea, Form, Select, ChoiceField,IntegerField,TextInput
from django.forms.widgets import TextInput
from django.forms.formsets import BaseFormSet

from designerGui.models import DomainOrder, SubstrateOrder
from designerGui.models import make_uuid

class NRPForm(ModelForm):
	description = CharField(widget=Textarea, required=False)
	class Meta:
		model = NRP
		exclude = ['uuid', 'monomers','owner','designed','construct','designerDomains']



def make_changed_boundary_nrp_form(nrp_uuid):
    class ChangedBoundaryNRPForm(Form):
        domain = ChoiceField(label="Domain", widget=Select())
        left_boundary = IntegerField(label="Left Boundary", widget=TextInput)
        right_boundary = IntegerField(label="Right Boundary", widget=TextInput)

        def __init__(self,  *args, **kwargs):
            super(Form, self).__init__(*args, **kwargs)
            self.nrp_uuid = nrp_uuid
            nrp = NRP.objects.get(uuid = nrp_uuid)
            self.fields['domain'].choices = [(x.pk, str(x)) for x in DomainOrder.objects.filter(nrp=nrp)]

        def save(self):
            # extract data from cleaned form
            chosen_dom_order = DomainOrder.objects.get(pk = self.cleaned_data['domain'])
            left_boundary    = self.cleaned_data['left_boundary']
            right_boundary   = self.cleaned_data['right_boundary']

            # start saving nrp variants
            nrp = NRP.objects.get(uuid = self.nrp_uuid)
            parent_dom_order_all = DomainOrder.objects.filter(nrp=nrp)
            parent_nrp_order_all = SubstrateOrder.objects.filter(nrp=nrp)
            nrp.pk = None
            nrp.name = str(nrp.name) + " Boundary variant"
            nrp.description = str(chosen_dom_order) + " Boundaries:" + str(left_boundary) + "-" + str(right_boundary)
            nrp.construct = None
            nrp.uuid = make_uuid()
            nrp.boundary_parent = NRP.objects.get(uuid=self.nrp_uuid)


            nrp.save()

            for nrp_order in parent_nrp_order_all:
                nrp_order.pk = None
                nrp_order.nrp = nrp
                nrp_order.save()

            for domain_order in parent_dom_order_all:
                if domain_order == chosen_dom_order:
                    domain_order.left_boundary = left_boundary
                    domain_order.right_boundary = right_boundary
                domain_order.pk = None
                domain_order.nrp = nrp
                domain_order.save()

            nrp.makeConstruct()
            NRP.objects.get(uuid=self.nrp_uuid).adjustConstruct(nrp)

    return ChangedBoundaryNRPForm


