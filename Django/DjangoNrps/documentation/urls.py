from django.conf.urls import patterns, url

from documentation.views import DocumentationView

urlpatterns = patterns('',
    url(r'^$', DocumentationView.as_view(), name='documentation'),
)