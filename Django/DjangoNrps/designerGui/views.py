from designerGui.models import Species
from databaseInput.models import Substrate
from databaseInput.forms import SubstrateFormSet
from django.views.generic import ListView, CreateView
from django.http import HttpResponse

import openbabel as ob

class SpeciesListView(ListView):
  template_name = 'designerGui/use_tool.html'
  model = Species
  
  def get_context_data(self, **kwargs):

        context = super(SpeciesListView, self).get_context_data(**kwargs)
        context['myFormSet'] = SubstrateFormSet()
        return context
  
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
    conv.ReadString(mol, str(Substrate.objects.get(pk=int(aminoacids[0])).structure))
    pattern.Match(mol)
    mollist = pattern.GetUMapList()[0]
    natom = mol.GetAtom(mollist[0])
    makeResidue(mol, 0, mollist)

    i = 1
    for aa in aminoacids[1:]:
        mol2 = ob.OBMol()
        conv.ReadString(mol2,str(Substrate.objects.get(pk=int(aa)).structure))
        pattern.Match(mol2)
        mollist = pattern.GetUMapList()[0]
        makeResidue(mol2, i, mollist)

        catom = mol2.GetAtom(mollist[2])
        oatom = mol2.GetAtom(mollist[4])

        molnatoms = mol.NumAtoms()
        mol += mol2

        builder.Connect(mol, natom.GetIdx(), molnatoms + catom.GetIdx())
        foatom = mol.GetAtom(molnatoms + oatom.GetIdx())
        fnatom = mol.GetAtom(molnatoms + mol2.GetAtom(mollist[0]).GetIdx())
        mol.DeleteHydrogens(foatom)
        mol.DeleteAtom(foatom)

        natom.SetImplicitValence(3)
        mol.DeleteHydrogens(natom)
        mol.AddHydrogens(natom)
        natom = fnatom
        i += 1

    builder.Build(mol)
    for res in ob.OBResidueIter(mol):
        for atom in ob.OBResidueAtomIter(res):
            for bond in ob.OBAtomBondIter(atom):
                data = ob.OBPairData()
                data.SetAttribute("color")
                data.SetValue(color)
                bond.CloneData(data)
    mol.DeleteHydrogens()
    svg = conv.WriteString(mol)
    return HttpResponse(svg, mimetype="image/svg+xml")

