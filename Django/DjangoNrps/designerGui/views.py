from designerGui.models import Species, NRP, SubstrateOrder
from databaseInput.models import Substrate, Modification, Domain, Type
from databaseInput.forms import SubstrateFormSet, ModificationsFormSet
from designerGui.forms import NRPForm, make_changed_boundary_nrp_form
from gibson.jsonresponses import JsonResponse, ERROR
from gibson.views import primer_download

from django.shortcuts import redirect
from django.views.generic import ListView, CreateView, TemplateView
from django.http import HttpResponse, HttpResponseNotFound, HttpResponseRedirect, QueryDict, Http404
from django.contrib.auth.decorators import login_required
from django.template import Context, loader, RequestContext
from django.core.urlresolvers import reverse
from django.views.decorators.cache import cache_page
from django.core.exceptions import ObjectDoesNotExist
from django.forms.formsets import formset_factory


import math
import json
import xml.etree.ElementTree as x
import openbabel as ob
import json
from uuid import uuid4

def toBool(x):
    return str(x.lower()) in ("yes", "true", "t", "1")

def get_nrp(user, uuid):
    if user.is_authenticated():
        owner = user
    else:
        owner = None
    try:
        nrp = NRP.objects.get(owner=owner, uuid=uuid)
    except ObjectDoesNotExist:
        return False
    else:
        return nrp

def makeConstruct(request, uuid):
    nrp = get_nrp(request.user, uuid)
    if not nrp:
        return HttpResponseNotFound()
    else:
        nrp.designed = False
        return JsonResponse({'taskId': nrp.designDomains.delay(toBool(request.POST['curatedonly'])).id})

def getConstruct(request, uuid):
    nrp = get_nrp(request.user, uuid)
    if not nrp:
        return HttpResponseNotFound()
    else:
        if not nrp.designed:
            return makeConstruct(request, uuid)
        elif nrp.construct is None or nrp.construct.fragments.count() == 0:
            con = nrp.makeConstruct()
        con = nrp.construct
        constructId = con.pk
        designTabLink = reverse('design_tab', kwargs= {'cid' : constructId})
        primerTabLink = reverse('primers', kwargs= {'cid' : constructId})
        domainSequenceTabLink = reverse('domainSequence', kwargs = {'uuid' : nrp.uuid})
        return JsonResponse({"constructId": constructId,
                            'designTabLink': designTabLink,
                            'primerTabLink': primerTabLink,
                            'domainSequenceTablLink': domainSequenceTabLink})

def nrpDesigner(request, uuid, sample=False):
    nrp = get_nrp(request.user, uuid)
    if not nrp and request.user.is_authenticated():
        return HttpResponseNotFound()
    if not nrp:
        if NRP.objects.filter(uuid=uuid).count() > 0:
            return nrpListView(request, sample)
        nrp = NRP(owner=None, uuid=uuid)
        if sample:
            nrp.save()
            nrp.makeSample()
    t = loader.get_template('gibson/designer.html')
    c = RequestContext(request, {'nrp': nrp})
    return HttpResponse(t.render(c))

@login_required
def peptide_add(request):
    if request.method == 'POST':
        form = NRPForm(request.POST, prefix='nrp')
        if form.is_valid():
            c = form.save(commit=False)
            c.owner = request.user
            c.save()
            return JsonResponse({'url': reverse('peptides') })
        t = loader.get_template('designerGui/peptideform.html')
        con = NRP.objects.all().filter(owner=request.user)
        c = RequestContext(request, {
            'NRPform':form,
        })
        return JsonResponse({'html': t.render(c),}, ERROR)
    else:
        return HttpResponseNotFound()

def peptide_delete(request, uuid):
    peptide = get_nrp(request.user, uuid)
    if peptide:
        peptide.fullDelete()
        if request.is_ajax():
            return JsonResponse('/peptides')
        return HttpResponseRedirect('/peptides')
    else:
        return HttpResponseNotFound()

def nrpListView(request, sample=False):
    if request.user.is_authenticated():
        t = loader.get_template('designerGui/peptides.html')
        c = RequestContext(request, {'NRPform': NRPForm(prefix='nrp'), 'title': 'NRPSDesigner', 'nrpList': NRP.objects.filter(owner=request.user)})
        return HttpResponse(t.render(c))
    else:
        uuid = str(uuid4())
        if not sample:
            return redirect('nrpDesigner', uuid)
        else:
            return redirect('sampleNrpDesigner', uuid)


def SpeciesListView(request, uuid):
    nrp = get_nrp(request.user, uuid)
    if not nrp and request.user.is_authenticated():
        return HttpResponseNotFound()
    if not nrp:
        nrp = NRP(owner=None, uuid=uuid)
    t = loader.get_template('designerGui/use_tool.html')
    substrateOrder = SubstrateOrder.objects.filter(nrp = nrp)
    modfs = Modification.objects.all()
    c = RequestContext(request, {'nrp': nrp, 'uuid': uuid, 'myFormSet': SubstrateFormSet(), 'modifications': modfs.values(), 'substrateOrder': substrateOrder, 'indigoidineTagged': nrp.indigoidineTagged, 'initialPic': nrp.getPeptideSequenceForStructView()})
    return HttpResponse(t.render(c))

def createLibrary(request, uuid):
    nrp = get_nrp(request.user, uuid)
    if not nrp:
        return HttpResponseNotFound()
    t = loader.get_template('designerGui/createlibrary.html')
    c = RequestContext(request)
    substrateOrder = SubstrateOrder.objects.filter(nrp=nrp)
    c['substrateOrder'] = substrateOrder
    c['indigoidineTagged'] = nrp.indigoidineTagged
    c['uuid'] = uuid
    initialPic = nrp.getPeptideSequenceForStructView()
    c['initialPic'] = initialPic
    c['scaffold'] = nrp.monomers.order_by('substrateOrder')
    return HttpResponse(t.render(c))

def processLibrary(request, uuid):
    if request.method == 'POST':
        nrp = get_nrp(request.user, uuid)
        if not nrp:
            return HttpResponseNotFound()
        data = json.loads(request.body)
        monomers = []
        haveUseAll = -1
        i = 0
        for m in data['as']:
            monomers.append(m[0])
            if len(m) > 1 and m[1] == True:
                haveUseAll = i
            i += 1
        if haveUseAll != -1:
            for m in data['as']:
                del m[1:]
            substrates = get_available_substrates(monomers, True, data['curatedonly'], haveUseAll)
            scaffold = data['as'][haveUseAll][0]
            for s in substrates:
                if s.pk != scaffold:
                    data['as'][haveUseAll].append(s.pk)
        return JsonResponse({'taskId': nrp.makeLibrary.delay(data['as'], data['curatedonly']).id})

def viewLibrary(request, uuid):
    parentnrp = get_nrp(request.user, uuid)
    if not parentnrp:
        return HttpResponseNotFound()
    childnrps = parentnrp.child.all()
    c = RequestContext(request)
    c['parentnrp'] = parentnrp
    c['childnrps'] = childnrps
    t = loader.get_template('designerGui/viewlibrary.html')
    return HttpResponse(t.render(c))

def downloadLibrary(request, uuid):
    parentnrp = get_nrp(request.user, uuid)
    if not parentnrp:
        return HttpResponseNotFound()
    cids = [child.construct.pk for child in parentnrp.child.filter(uuid__in=request.POST.getlist('id'))]
    return primer_download(request, parentnrp.construct.pk, *cids)

def create_boundary_library(request, uuid):
    nrp = get_nrp(request.user, uuid)
    if not nrp:
        return HttpResponseNotFound()

    ChangedBoundaryNRPForm = make_changed_boundary_nrp_form(uuid)
    changedBoundaryNRPFormset = formset_factory(ChangedBoundaryNRPForm, extra=1)

    if request.method == "POST":
        boundary_formset = changedBoundaryNRPFormset(request.POST)
        if boundary_formset.is_valid():
            for boundary_form in boundary_formset:
                if boundary_form.is_valid():
                    boundary_form.save()
            return HttpResponseRedirect(reverse("viewBoundaryLibrary", kwargs={'uuid': uuid}))
        else:
            t = loader.get_template('designerGui/createBoundaryLibrary.html')
            c = RequestContext(request,{
                    'formset': boundary_formset
                })

        return HttpResponse(t.render(c))

    else:
        #boundary_form = ChangedBoundaryNRPForm(nrp)
        t = loader.get_template('designerGui/createBoundaryLibrary.html')
        c = RequestContext(request,{
                    'formset': changedBoundaryNRPFormset
                })

        return HttpResponse(t.render(c))

def view_boundary_library(request, uuid):
    parentnrp = get_nrp(request.user, uuid)
    if not parentnrp:
        return HttpResponseNotFound()
    childnrps = parentnrp.boundary_child.all()
    c = RequestContext(request)
    c['parentnrp'] = parentnrp
    c['childnrps'] = childnrps
    c['libraryclassification'] = 'boundary'
    t = loader.get_template('designerGui/viewlibrary.html')
    return HttpResponse(t.render(c))

def get_available_monomers(request):
    if request.method == 'POST' and "monomer[]" in request.POST:
        monomers = request.POST.getlist("monomer[]")
        if 'selected' in request.POST:
            selected = int(request.POST['selected'])
        else:
            selected = None
        aas = get_available_substrates(monomers, toBool(request.POST['current']), toBool(request.POST['curatedonly']), selected)
    else:
        aas = filter(lambda x: x.can_be_added(), Substrate.objects.exclude(user__username='sbspks'))
    json = {}
    minid = float("Inf")
    for aa in aas:
        if aa.parent is None:
            name = aa.name
            if aa.name[0:2].upper() == 'L-' or aa.name[0:2].upper() == 'D-':
                name = aa.name[2:]
            if aa.pk in json:
                key = aa.pk
            elif aa.enantiomer is not None and aa.enantiomer.pk in json:
                key = aa.enantiomer.pk
            else:
                key = None
            if aa.enantiomer is None:
                chirality = 'N'
            else:
                chirality = aa.chirality
            if key is not None:
                if key < minid:
                    minid = key
                json[key][chirality.lower() + "id"] = aa.pk
                json[key][chirality+'Children'] = [{"text": c.name, "id": c.pk} for c in aa.child.all()]
                #names[name]['name'] = name
            else:
                json[aa.pk] = {"id": aa.pk, chirality.lower() + "id": aa.pk, 'text': name, aa.chirality+'Children': [{"text": c.name, "id": c.pk} for c in aa.child.all()]}

    jsonlist = json.values()
    jsonlist.sort(lambda x,y: cmp(x['text'], y['text']))
    return JsonResponse({"monomers": json, "monomerslist": jsonlist})

def get_available_substrates(monomers, current, curatedonly=True, selected=None):
    if selected is None:
        selected = len(monomers) - 1
    if current:
        if selected > 0:
            chirality = Substrate.objects.get(pk=monomers[selected - 1]).chirality
        else:
            chirality = None
    else:
        chirality = Substrate.objects.get(pk=monomers[-1]).chirality
    substrates = Substrate.objects.exclude(user__username='sbspks')
    if not current or selected == len(monomers) - 1:
        aas = filter(lambda x: x.can_be_added(chirality, curatedonly), substrates)
    else:
        following = Substrate.objects.get(pk=monomers[selected + 1])
        aas = filter(lambda x: x.can_be_added(chirality, curatedonly) and following.can_be_added(x.chirality, curatedonly), substrates)
    return aas

@login_required
def submit_nrp(request):
    nrpxml = x.Element('nrp')
    for monomer in request.POST.getlist("as[]"):
        monomerel = x.SubElement(nrpxml, 'monomer')
        monomerid = x.SubElement(monomerel, 'id')
        monomerid.text = monomer
    return HttpResponse(x.tostring(nrpxml, "utf8"), mimetype="text/xml")

@cache_page(365*24*60*60*15)
def make_structure(request):
    def makeResidue(mol, idx, aaatoms):
        res = mol.NewResidue()
        res.SetNum(idx)
        for atom in ob.OBMolAtomIter(mol):
            if atom.GetIdx() not in aaatoms:
                res.AddAtom(atom)

    aminoacids = request.GET.getlist("as")
    color = request.GET.get("rescol", "#3EC1CD")

    mol = ob.OBMol()
    conv = ob.OBConversion()
    pattern = ob.OBSmartsPattern()
    pattern.Init("[NX3][$([CX4H1]([*])),$([CX4H2])][CX3](=[OX1])[OX2]")
    builder = ob.OBBuilder()
    conv.SetInAndOutFormats("sdf", "svg")
    conv.AddOption("d", ob.OBConversion.OUTOPTIONS)
    conv.AddOption("b", ob.OBConversion.OUTOPTIONS, "none")
    conv.ReadString(mol, str(Substrate.objects.get(pk=int(aminoacids[0])).structure))
    pattern.Match(mol)
    mollist = pattern.GetUMapList()[0]
    oatom = mol.GetAtom(mollist[4])
    catom = mol.GetAtom(mollist[2])
    firstnatom = mol.GetAtom(mollist[0])
    makeResidue(mol, 0, mollist)

    i = 1
    for aa in aminoacids[1:]:
        mol2 = ob.OBMol()
        conv.ReadString(mol2, str(Substrate.objects.get(pk=int(aa)).structure))
        pattern.Match(mol2)
        mollist = pattern.GetUMapList()[0]
        makeResidue(mol2, i, mollist)

        molnatoms = mol.NumAtoms()
        mol += mol2

        natom = mol.GetAtom(molnatoms + mol2.GetAtom(mollist[0]).GetIdx())

        builder.Connect(mol, catom.GetIdx(), natom.GetIdx())

        foatom = mol.GetAtom(molnatoms + mol2.GetAtom(mollist[4]).GetIdx())
        catom = mol.GetAtom(molnatoms + mol2.GetAtom(mollist[2]).GetIdx())
        mol.DeleteHydrogens(oatom)
        mol.DeleteAtom(oatom)

        natom.SetImplicitValence(3)
        mol.DeleteHydrogens(natom)
        mol.AddHydrogens(natom)
        oatom = foatom

        i += 1

    nidx = firstnatom.GetIdx()
    oidx =  oatom.GetIdx()
    builder.Build(mol)
    natom = mol.GetAtom(nidx)
    oatom = mol.GetAtom(oidx)
    for res in ob.OBResidueIter(mol):
        for atom in ob.OBResidueAtomIter(res):
            for bond in ob.OBAtomBondIter(atom):
                data = ob.OBPairData()
                data.SetAttribute("color")
                data.SetValue(color)
                bond.CloneData(data)
    mol.DeleteHydrogens()
    gen2d = ob.OBOp.FindType("gen2d")
    gen2d.Do(mol)

    opp = oatom.GetY() - natom.GetY()
    adj = oatom.GetX() - natom.GetX()
    angle = abs(math.atan(opp / adj))
    if opp > 0 and adj > 0:
        pass
    elif opp > 0 and adj < 0:
        angle = math.pi - angle
    elif opp < 0 and adj < 0:
        angle = math.pi + angle
    elif opp < 0 and adj > 0:
        angle = 2 * math.pi - angle
    angle = -angle
    mol.Rotate(ob.double_array([math.cos(angle), -math.sin(angle), 0,
                                math.sin(angle), math.cos(angle),  0,
                                0,               0,                1]))
    svg = conv.WriteString(mol)
    # need to get rid of square aspect ratio
    delstart = svg.find("width")
    delend = svg.find("svg", delstart)
    delend = svg.find("viewBox", delend)
    svgend = svg.rfind("</g>")
    svg = svg[0:delstart] + svg[delend:svgend]
    return HttpResponse(svg, mimetype="image/svg+xml")

def domainSequenceView(request, uuid):
    nrp = get_nrp(request.user, uuid)
    if not nrp:
        return HttpResponseNotFound()
    t = loader.get_template('designerGui/domainSequence.html')
    context = RequestContext(request, {'pfamGraphicJson': nrp.generatePfamGraphicJson()})
    return HttpResponse(t.render(context))
