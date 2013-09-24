from django.conf.urls import patterns, include, url
from designerGui.views import SpeciesListView, make_structure, submit_nrp, NRPListView, peptide_add, peptide_delete, nrpDesigner, makeConstruct, DomainSequenceView
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse, HttpResponseRedirect

from designerGui.api import saveNrpMonomers

urlpatterns = patterns('',
					url(r'^$', login_required(NRPListView.as_view()),name="peptides"),
					url(r'^add$', peptide_add, name="peptide_add"),
					url(r'^(?P<cid>\d+)/delete/$', peptide_delete, name = "peptide_delete"),
					url(r'^(?P<pid>\d+)/$', nrpDesigner, name = "nrpDesigner"),
					url(r'^nrpselection/(?P<pid>\d+)$', SpeciesListView.as_view(), name="guiTool"),
					url(r'^makeConstruct/(?P<pid>\d+)$', makeConstruct, name="makeConstruct"),
					url(r'^structure/', make_structure, name="Structure"),
					url(r'^submit/', submit_nrp, name="submitNRP"),
					url(r'^api/(?P<pid>\d+)', saveNrpMonomers, name = "saveNrpMonomers" ),
					url(r'^(?P<pid>\d+)/domainSequence/', DomainSequenceView.as_view(), name="domainSequence")
					)
