from django.http import HttpResponse, HttpResponseRedirect, HttpResponseServerError, HttpResponseBadRequest
from django.conf import settings
from django.core.exceptions import ObjectDoesNotExist
from databaseInput.models import Origin, Type, Domain, Cds, ORIGIN_TYPES
from designerGui.models import DomainOrder
from gibson.jsonresponses import JsonResponse, ERROR

def get_origins(request):
    jsondata = []
    for otype in ORIGIN_TYPES:
        origins = Origin.objects.filter(sourceType=otype[0])
        if len(origins) > 0:
            data = [{'id': o.pk, 'text': str(o)} for o in origins]
            data.sort(lambda x,y: cmp(x['text'], y['text']))
            jsondata.append({'text': origins[0].get_sourceType_display(), 'children': data})
    return JsonResponse(jsondata)

def get_cds(request):
    if request.method == "GET":
        origin_id = request.GET['originId']
        if origin_id == "":
            cds = Cds.objects.all()
        else:
            cds = Cds.objects.filter(origin=origin_id)
        jsondata = [{'id': c.pk, 'text': c.geneName} for c in cds]
        jsondata.sort(lambda x,y: cmp(x['text'], y['text']))
        return JsonResponse(jsondata)

def get_domain_types(request):
    domain_types = Type.objects.all()
    ret = [{'id':x.pk, 'text':str(x)} for x in domain_types]
    return JsonResponse(ret)

def get_info(request):
    if request.method == "POST":
        jsondata = []
        try:
            if 'cdsId' in request.POST and ('domainId' not in request.POST or len(request.POST['domainId']) == 0):
                domain = None
                cds = Cds.objects.get(pk=request.POST['cdsId'])
            elif 'domainId' in request.POST and len(request.POST['domainId']) > 0:
                domain = Domain.objects.get(pk=request.POST['domainId'])
                cds = domain.cds
            else:
                return HttpResponseNotFound()
        except ObjectDoesNotExist:
            return HttpResponseNotFound()
        ori = cds.origin
        species = ori.species
        while species is None and ori.parent is not None:
            ori = ori.parent
            species = ori.species

        jsondata.append(('Species', species))
        jsondata.append(('Gene', cds.geneName))
        jsondata.append(('Product', cds.product.name))
        if cds.origin.sourceType == "Species":
            jsondata.append(('Taxonomy ID', cds.origin.source))
        elif cds.origin.sourceType == "Biobrick":
            jsondata.append(('BioBrick', cds.origin.source))
        else:
            jsondata.append(('Source', cds.origin.source))
        if domain is None:
            jsondata.append(('Description', cds.description))
            user = cds.user
        else:
            jsondata.append(('Domain', domain.domainType.name))
            jsondata.append(('Module', domain.module))
            if domain.substrateSpecificity is not None:
                substrates = [s.name for s in domain.substrateSpecificity.all()]
                jsondata.append(('Substrate', ", ".join(substrates)))
            if domain.chirality != 'N':
                jsondata.append(('Chirality', domain.chirality))
            jsondata.append(('Description', domain.description))
            user = domain.user
        if user is not None and user.groups.filter(name=settings.CURATION_GROUP).count() > 0:
            curated = "Yes"
        else:
            curated = "No"
        jsondata.append(('Curated', curated))
        return JsonResponse(jsondata)
    else:
        return HttpResponseNotFound()


def get_domains(request):
    if request.method == "GET":
        if 'originId' in request.GET:
            origin_id = request.GET['originId']
        else:
            origin_id = ""
        if 'typeId' in request.GET:
            type_id = request.GET['typeId']
        else:
            type_id = ""
        if 'cdsId' in request.GET:
            cds_id = request.GET['cdsId']
        else:
            cds_id = ""
        if len(origin_id) == 0 and len(type_id) == 0 and len(cds_id) == 0:
            domains = Domain.objects.all()
        else:
            filters = {}
            if len(origin_id) > 0:
                filters['cds__origin__pk'] = origin_id
            if len(type_id) > 0:
                filters['domainType__pk'] = type_id
            if len(cds_id) > 0:
                filters['cds__pk'] = cds_id
            domains = Domain.objects.filter(**filters)
        ret = [{'id':x.pk, 'text':x.short_name()} for x in domains]
        return JsonResponse(ret)

def get_biojs_sequence(request):
    if request.method == "GET":
        try:
            jsondata = {}
            if 'cdsId' in request.GET:
                cds_id = request.GET['cdsId']
                cds = Cds.objects.get(pk=cds_id)
                jsondata['biojs'] = cds.get_biojs_entry()
            if 'domainId' in request.GET:
                domain_id = request.GET['domainId']
                domain = Domain.objects.get(pk = domain_id)
                if 'domainOnly' not in request.GET or request.GET['domainOnly'] == "false":
                    jsondata['biojs'] = domain.cds.get_biojs_entry()
                jsondata['highlight'] = domain.get_biojs_highlight()
            return JsonResponse(jsondata)
        except ObjectDoesNotExist:
            return HttpResponseNotFound()
    else:
        return HttpResponseNotFound()

def get_domain_order_biojs_sequence(request):
    if request.method == "GET" and 'domainOrderId' in request.GET:
        try:
            jsondata = {}
            domain_order_id = request.GET['domainOrderId']
            domain_order = DomainOrder.objects.get(pk = domain_order_id)
            domain = domain_order.domain
            jsondata['biojs'] = domain.cds.get_biojs_entry(protein=True)
            jsondata['highlight'] = domain.get_biojs_highlight(protein=True)
            return JsonResponse(jsondata)
        except ObjectDoesNotExist:
            return HttpResponseNotFound()
    else:
        return HttpResponseNotFound()
