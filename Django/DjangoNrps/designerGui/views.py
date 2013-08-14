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
    aminoacids = request.GET.getlist("as")

    mol = ob.OBMol()
    conv = ob.OBConversion()
    pattern = ob.OBSmartsPattern()
    pattern.Init("[NX3][$([CX4H1]([*])),$([CX4H2])][CX3](=[OX1])[OX2]")
    builder = ob.OBBuilder()
    conv.SetInAndOutFormats("sdf", "svg")
    conv.ReadString(mol, str(Substrate.objects.get(pk=int(aminoacids[0])).structure))
    pattern.Match(mol)
    mollist = pattern.GetUMapList()
    natom = mol.GetAtom(mollist[0][0])

    for aa in aminoacids[1:]:
        mol2 = ob.OBMol()
        conv.ReadString(mol2,str(Substrate.objects.get(pk=int(aa)).structure))
        pattern.Match(mol2)
        mollist = pattern.GetUMapList()

        catom = mol2.GetAtom(mollist[0][2])
        oatom = mol2.GetAtom(mollist[0][4])

        molnatoms = mol.NumAtoms()
        mol += mol2

        builder.Connect(mol, natom.GetIdx(), molnatoms + catom.GetIdx())
        foatom = mol.GetAtom(molnatoms + oatom.GetIdx())
        fnatom = mol.GetAtom(molnatoms + mol2.GetAtom(mollist[0][0]).GetIdx())
        mol.DeleteHydrogens(foatom)
        mol.DeleteAtom(foatom)
        #mol.DeleteHydrogens()
        #mol.AddHydrogens()

        natom.SetImplicitValence(3)
        mol.DeleteHydrogens(natom)
        mol.AddHydrogens(natom)
        natom = fnatom

    builder.Build(mol)
    svg = conv.WriteString(mol)
    return HttpResponse(svg, mimetype="image/svg+xml")

