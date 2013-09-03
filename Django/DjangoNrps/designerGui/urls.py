from django.conf.urls import patterns, include, url
from designerGui.views import SpeciesListView, make_structure, submit_nrp, NRPListView, peptide_add, peptide_delete
from django.contrib.auth.decorators import login_required

urlpatterns = patterns('',
					url(r'^peptides/$', login_required(NRPListView.as_view()),name="peptides"),
					url(r'^peptides/add$', peptide_add, name="peptide_add"),
					url(r'^peptides/(?P<cid>\d+)/delete/$', peptide_delete, name = "peptide_delete"),
					url(r'^$', SpeciesListView.as_view(), name="guiTool"),
					url(r'^structure/', make_structure, name="Structure"),
					url(r'^submit/', submit_nrp, name="submitNRP"),
					)
