from django.conf.urls import patterns, url

from documentation.views import DocumentationView, VideoTutorialView

urlpatterns = patterns('',
	url(r'^video/$', VideoTutorialView.as_view(), name='video_tutorial'),
	url(r'^$', DocumentationView.as_view(), name='documentation'),
)