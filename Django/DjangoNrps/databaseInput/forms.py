from databaseInput.models import Origin, Cds, Domain
from django.forms.models import inlineformset_factory


CdsFormSet = inlineformset_factory(Origin, Cds, extra=1)
DomainFormSet = inlineformset_factory(Cds, Domain, extra=1)
