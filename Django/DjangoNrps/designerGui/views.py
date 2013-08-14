from designerGui.models import Species
from databaseInput.models import Substrate
from databaseInput.forms import SubstrateFormSet
from django.views.generic import ListView, CreateView


class SpeciesListView(ListView):
  template_name = 'designerGui/use_tool.html'
  model = Species
  
  def get_context_data(self, **kwargs):

        context = super(SpeciesListView, self).get_context_data(**kwargs)
        context['myFormSet'] = SubstrateFormSet()
        return context
  
