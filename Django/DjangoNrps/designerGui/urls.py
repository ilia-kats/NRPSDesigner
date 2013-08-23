from django.conf.urls import patterns, include, url
from designerGui.views import SpeciesListView, make_structure, submit_nrp

urlpatterns = patterns('',
                       url(r'^$', SpeciesListView.as_view(), name="guiTool"),
                       url(r'^structure/', make_structure, name="Structure"),
                       url(r'^submit/', submit_nrp, name="submitNRP"),
                       )
