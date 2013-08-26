from designerGui.models import Species
from databaseInput.models import Substrate, Modification
from databaseInput.forms import SubstrateFormSet, ModificationsFormSet
from django.views.generic import ListView, CreateView
from django.http import HttpResponse

import math
import json
import xml.etree.ElementTree as x
import openbabel as ob

class SpeciesListView(ListView):
  template_name = 'designerGui/use_tool.html'
  model = Species
  
  def get_context_data(self, **kwargs):

        context = super(SpeciesListView, self).get_context_data(**kwargs)
        context['myFormSet'] = SubstrateFormSet()
        
        modfs = Modification.objects.all()
        context['modifications'] = modfs.values()
        #context['modifications'].sort(lambda x,y: cmp(x['name'], y['name']))
        
        
        aas = Substrate.objects.all()
        realAas = []
        for aa in aas:
            if not hasattr(aa.parent, 'name'):
                realAas.append(aa)
        names = {}
        for aa in realAas:
            name = aa.name
            if aa.name[0:2].upper() == 'L-' or aa.name[0:2].upper() == 'D-':
                name = aa.name[2:]
            if name in names:
                names[name][aa.chirality] = aa.pk
                names[name][aa.chirality+'Children'] = aa.child.all()
                #names[name]['name'] = name
            else:
                tmp = aa.child.all()
                names[name] = {aa.chirality: aa.pk, 'name': name, aa.chirality+'Children': aa.child.all()}
        context['substrates'] = names.values()
        context['substrates'].sort(lambda x,y: cmp(x['name'], y['name']))
        return context

def submit_nrp(request):
    nrpxml = x.Element('nrp')
    for monomer in request.POST.getlist("as[]"):
        monomerel = x.SubElement(nrpxml, 'monomer')
        monomerid = x.SubElement(monomerel, 'id')
        monomerid.text = monomer
    return HttpResponse(x.tostring(nrpxml, "utf8"), mimetype="text/xml")
  
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
