from databaseInput.models import Origin
from django.views.generic import ListView


class OriginListView(ListView):
  template_name = 'designerGui/use_tool.html'
  model = Origin