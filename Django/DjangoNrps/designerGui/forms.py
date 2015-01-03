from designerGui.models import NRP
from django.forms.models import inlineformset_factory, formset_factory
from django.forms import ModelForm, CharField, Textarea, Form, Select, ChoiceField,IntegerField,TextInput, HiddenInput
from django.forms import ValidationError
from django.forms.widgets import TextInput
from django.forms.formsets import BaseFormSet

from designerGui.models import DomainOrder, SubstrateOrder
from designerGui.models import make_uuid
from databaseInput.models import prot_to_dna_coords, dna_to_prot_coords

from ast import literal_eval

class NRPForm(ModelForm):
	description = CharField(widget=Textarea, required=False)
	class Meta:
		model = NRP
		exclude = ['uuid', 'monomers','owner','designed','construct','designerDomains']



def make_changed_boundary_nrp_form(nrp_uuid):
    class ChangedBoundaryNRPForm(Form):
        parent = CharField(required=False, widget=HiddenInput())
        linker = ChoiceField(label="Linker", widget=Select())
        left_boundary = IntegerField(label="Left Boundary", widget=TextInput)
        right_boundary = IntegerField(label="Right Boundary", widget=TextInput)

        def __init__(self,  *args, **kwargs):
            super(Form, self).__init__(*args, **kwargs)
            self.nrp_uuid = nrp_uuid
            nrp = NRP.objects.get(uuid = nrp_uuid)
            domains = DomainOrder.objects.filter(nrp=nrp)
            self.fields['linker'].choices = [((x.pk, y.pk), "%s - %s" % (str(x), str(y))) for x,y in zip(domains, domains[1:])]

        def save(self, parent=None):
            # extract data from cleaned form
            linker = literal_eval(self.cleaned_data.get("linker"))
            left_domain = DomainOrder.objects.get(pk = linker[0])
            right_domain = DomainOrder.objects.get(pk=linker[1])
            left_boundary    = self.cleaned_data['left_boundary']
            right_boundary   = self.cleaned_data['right_boundary']


            # start saving nrp variants
            if parent is None:
                nrp = NRP.objects.get(uuid = self.nrp_uuid)
            else:
                nrp = parent
            parent_dom_order_all = DomainOrder.objects.filter(nrp=nrp)
            parent_nrp_order_all = SubstrateOrder.objects.filter(nrp=nrp)
            nrp.pk = None
            nrp.name = str(nrp.name)
            desc = str("%s - %s" % (str(left_domain), str(right_domain))) + " Boundaries:" + str(left_boundary) + "-" + str(right_boundary) + str("aa")
            if nrp.name.endswith(" Boundary variant"):
                nrp.description += "%s " % desc
            else:
                nrp.name += " Boundary variant"
                nrp.description = desc
            nrp.description
            nrp.construct = None
            nrp.uuid = make_uuid()
            nrp.boundary_parent = NRP.objects.get(uuid=self.nrp_uuid)


            nrp.save()

            for nrp_order in parent_nrp_order_all:
                nrp_order.pk = None
                nrp_order.nrp = nrp
                nrp_order.save()

            left_boundary, right_boundary = prot_to_dna_coords(left_boundary,right_boundary)

            for domain_order in parent_dom_order_all:
                if domain_order == left_domain:
                    domain_order.right_boundary = left_boundary
                elif domain_order == right_domain:
                    domain_order.left_boundary = right_boundary
                domain_order.pk = None
                domain_order.nrp = nrp
                domain_order.save()

            nrp.makeConstruct()
            NRP.objects.get(uuid=self.nrp_uuid).adjustConstruct(nrp)
            return nrp

        def clean(self):
            cleaned_data = super(ChangedBoundaryNRPForm, self).clean()
            left_boundary = cleaned_data.get("left_boundary")
            right_boundary = cleaned_data.get("right_boundary")
            linker = literal_eval(cleaned_data.get("linker"))
            left_domain = DomainOrder.objects.get(pk = linker[0])
            right_domain = DomainOrder.objects.get(pk=linker[1])

            left_min,right_max = dna_to_prot_coords(1, len(right_domain.domain.cds.dnaSequence))
            if  left_boundary < left_min:
                raise ValidationError("Left boundary has to be greater or equal to " + str(left_min)+".")
            if right_boundary > right_max:
                raise ValidationError("Right boundary exceeds length of sequence (" + str(right_max) +")")
            return cleaned_data

    return ChangedBoundaryNRPForm


