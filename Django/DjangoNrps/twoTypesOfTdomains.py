from databaseInput.models import Domain, Type

tDomains = Domain.objects.filter(domainType__name="T")

#lets get our juicy new domain types
tStd = Type.objects.get(name="Tstd")
tEp = Type.objects.get(name="T_ep")

# let's assume that Epimerization domains and preceding T-domain
# are always on the same coding sequence
for tDomain in tDomains:
	cds_domains = tDomain.cds.get_ordered_domain_list()
	possibleEDomainIndex = cds_domains.index(tDomain) + 1
	if possibleEDomainIndex < len(cds_domains) and cds_domains[possibleEDomainIndex].domainType.name == 'E':
		tDomain.domainType = tEp
	else:
		tDomain.domainType = tStd
	tDomain.save()

