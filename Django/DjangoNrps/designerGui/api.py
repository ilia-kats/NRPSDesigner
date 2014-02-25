from django.contrib.auth.decorators import login_required
from django.http import HttpResponse, HttpResponseNotFound
from designerGui.models import NRP, SubstrateOrder
from databaseInput.models import Substrate
from designerGui.views import toBool, get_nrp

def saveNrpMonomers(request, uuid):
    if request.method == 'POST':
        nrp = get_nrp(request.user, uuid)
        if not nrp and request.user.is_authenticated():
            return HttpResponseNotFound()
        if not nrp:
            nrp = NRP.objects.create(uuid=uuid, owner=None)
        else:
            # delete other stuff pointing
            nrp.delete_dependencies()
        # now add the new list
        for count, monomerId in enumerate(request.POST.getlist('as[]')):
            monomer = Substrate.objects.get(pk=int(monomerId))
            so = SubstrateOrder.objects.create(nrp= nrp, substrate = monomer, order = count)
        if 'indtag' in request.POST:
            nrp.indigoidineTagged = toBool(request.POST['indtag'])
        else:
            nrp.indigoidineTagged = False
        nrp.save()
        return HttpResponse()