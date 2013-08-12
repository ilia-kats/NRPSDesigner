from designerGui.models import Species
from django.views.generic import ListView


class SpeciesListView(ListView):
  template_name = 'designerGui/use_tool.html'
  model = Species