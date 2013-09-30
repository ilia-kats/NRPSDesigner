from django.contrib.auth.decorators import login_required
from django.http import HttpResponse
from designerGui.models import NRP, SubstrateOrder
from databaseInput.models import Substrate
from designerGui.views import toBool

@login_required
def saveNrpMonomers(request,pid):
    if request.method == 'POST':
        nrp = NRP.objects.get(owner = request.user, pk = pid)
        # delete other stuff pointing
        nrp.delete_dependencies()
        # now add the new list
        for count, monomerId in enumerate(request.POST.getlist('as[]')):
            monomer = Substrate.objects.get(pk=int(monomerId))
            so = SubstrateOrder.objects.create(nrp= nrp, substrate = monomer, order = count)
        nrp.indigoidineTagged = toBool(request.POST['indtag'])
        nrp.save()
        return HttpResponse()