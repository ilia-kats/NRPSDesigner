from django.conf.urls.defaults import patterns, include, url
from django.views.generic import TemplateView
from django.conf import settings

primerpatterns = patterns('gibson.views',
	url(r'^$', 'primers', name='primers'),
	(r'^download$', 'primer_download'),
	(r'^pdf$', 'pdf'),
	(r'^save$', 'primer_save'),
	(r'^reset$', 'primer_reset'),
	(r'^csv$', 'csv_primers'),
	url(r'^(?P<pid>\d+)/boxplot/$', TemplateView.as_view(template_name = 'gibson/boxplot.html'), name="primer_boxplot"),
	(r'^(?P<pid>\d+)/setLength/$', 'primer_length'),
	(r'^(?P<pid>\d+)/[\w\-_]+/$', 'primer'),
	(r'^(?P<pid>\d+)/$', 'load_primer'),
	(r'^(?P<pid>\d+)/[\w\-_]+/offset$', 'primer_offset'),
	(r'^\d+/[\w\-_]+/reset$', 'primer_reset'),
)

constructpatterns = patterns('gibson.designer',
	#(r'^$', 'designer'),
	url(r'^design/$', 'design_tab', name='design_tab'),
	(r'^settings/$', 'construct_settings'),
)

apipatterns = patterns('gibson.api',
	(r'^saveSettings/$', 'save_settings'),
	(r'^saveMeta/$', 'save_meta'),
	(r'^getInfo/$', 'get_info'),
	(r'^addFragment/$', 'fragment_add'),
	(r'^rmFragment/$', 'fragment_remove'),
	(r'^saveOrder/$', 'save_order'),
)

urlpatterns = patterns('gibson.views',
	(r'^$', 'constructs'),
	(r'^add$', 'construct_add'),
	(r'^(?P<cid>\d+)/\w+\.gb', 'download'),
	(r'^(?P<cid>\d+)/delete/$', 'construct_delete'),
	(r'^(?P<cid>\d+)/fragments/$', 'construct_fragment'),
	url(r'^(?P<cid>\d+)/process/$', 'process', name = 'process'),
	(r'^(?P<cid>\d+)/clipping/(?P<cfid>\d+)/$', 'fragment_clipping'),
	(r'^(?P<cid>\d+)/clipping/(?P<cfid>\d+)/apply/$', 'apply_clipping'),
	(r'^(?P<cid>\d+)/summary/$', 'summary'),
	(r'^(?P<cid>\d+)/primers/', include(primerpatterns)),
	(r'^api/(?P<cid>\d+)/', include(apipatterns)),
	(r'^(?P<cid>\d+)/', include(constructpatterns)),
	(r'^walkthrough/$',TemplateView.as_view(template_name = 'gibson/walkthrough.html')),
)

