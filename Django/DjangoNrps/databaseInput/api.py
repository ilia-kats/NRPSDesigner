from django.http import HttpResponse, HttpResponseRedirect, HttpResponseServerError, HttpResponseBadRequest
from databaseInput.models import Origin, Type, Domain

from gibson.jsonresponses import JsonResponse, ERROR

def get_origins(request):
	origins = Origin.objects.all()
	ret = [{'id':x.pk, 'text':str(x)} for x in origins]
	return JsonResponse(ret)

def get_domain_types(request):
	domain_types = Type.objects.all()
	ret = [{'id':x.pk, 'text':str(x)} for x in domain_types]
	return JsonResponse(ret)

def get_domains(request):
	if request.method == "GET":
		origin_id = request.GET['originId']
		type_id   = request.GET['typeId']
		if origin_id == "" and type_id == "":
			domains = Domain.objects.all()
		elif type_id=="":
			domains = Domain.objects.filter(cds__origin__pk = origin_id)
		elif origin_id =="":
			domains = Domain.objects.filter(domainType__pk = type_id)
		else:
			domains = Domain.objects.filter(cds__origin__pk = origin_id, domainType__pk = type_id)
		ret = [{'id':x.pk, 'text':x.short_name()} for x in domains]
		return JsonResponse(ret)
