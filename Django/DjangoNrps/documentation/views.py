from django.views.generic import TemplateView

class DocumentationView(TemplateView):
    template_name = "documentation/documentation.html"

class VideoTutorialView(TemplateView):
    template_name = "documentation/video_tutorial.html"

