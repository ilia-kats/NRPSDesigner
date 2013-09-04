from django.conf.urls import patterns, include, url
from designerGui.views import SpeciesListView, make_structure, submit_nrp, NRPListView, peptide_add, peptide_delete, nrpDesign
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse, HttpResponseRedirect

from designerGui.api import saveNrpMonomers

urlpatterns = patterns('',
					url(r'^$', login_required(NRPListView.as_view()),name="peptides"),
					url(r'^add$', peptide_add, name="peptide_add"),
					url(r'^(?P<cid>\d+)/delete/$', peptide_delete, name = "peptide_delete"),
					#url(r'^(?P<cid>\d+)$', nrpDesign, name = "nrpDesign"),
					url(r'^(?P<pid>\d+)$', SpeciesListView.as_view(), name="guiTool"),
					url(r'^structure/', make_structure, name="Structure"),
					url(r'^submit/', submit_nrp, name="submitNRP"),
					url(r'^api/(?P<pid>\d+)', saveNrpMonomers, name = "saveNrpMonomers" )
					)
