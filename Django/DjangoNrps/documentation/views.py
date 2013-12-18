from django.views.generic import TemplateView

class DocumentationView(TemplateView):
    template_name = "documentation/documentation.html"
