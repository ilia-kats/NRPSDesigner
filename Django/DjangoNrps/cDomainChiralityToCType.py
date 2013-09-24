from databaseInput.models import Domain, Type

cDomains = Domain.objects.filter(domainType__name = 'C')
domTypes = {'L': Type.objects.get(name = 'C_L'),
	'D': Type.objects.get(name = 'C_D')}
for domain in cDomains:
	domain.domainType = domTypes[domain.chirality]
	domain.save()