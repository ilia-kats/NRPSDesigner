from databaseInput.models import Domain, Type

sbsDomains = Domain.objects.filter(user__username = "sbspks")

for domain in sbsDomains:
	domain.pfamLinkerStart = 3*domain.pfamLinkerStart-2
	domain.pfamLinkerStop = 3*domain.pfamLinkerStop-2
	domain.definedLinkerStart = 3*domain.definedLinkerStart-2
	domain.definedLinkerStop = 3*domain.definedLinkerStop-2
	domain.pfamStart = 3*domain.pfamStart - 2
	domain.pfamStop = 3*domain.pfamStop -2
	domain.definedStart = 3* domain.definedStart - 2
	domain.definedStop =  3*domain.definedStop -2
	domain.save()