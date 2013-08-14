from designerGui.models import Species
from databaseInput.models import Substrate
from django.views.generic import ListView


class SpeciesListView(ListView):
  template_name = 'designerGui/use_tool.html'
  model = Species
  
  def get_context_data(self, **kwargs):

        context = super(SpeciesListView, self).get_context_data(**kwargs)
        context['substrates'] = Substrate.objects.all()
        return context
  
