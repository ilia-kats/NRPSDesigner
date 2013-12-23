from django.conf.urls import patterns, include, url
from designerGui.views import SpeciesListView, make_structure, submit_nrp, nrpListView, peptide_add, peptide_delete, nrpDesigner, getConstruct, makeConstruct, domainSequenceView, get_available_monomers, createLibrary, processLibrary, viewLibrary, downloadLibrary
from django.contrib.auth.decorators import login_required
from django.http import HttpResponse, HttpResponseRedirect

from designerGui.api import saveNrpMonomers

uuidregex = r'[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}'

urlpatterns = patterns('',
					url(r'^$', nrpListView,name="peptides"),
					url(r'^add$', peptide_add, name="peptide_add"),
					url(r'^(?P<uuid>' + uuidregex + ')/delete/$', peptide_delete, name = "peptide_delete"),
					url(r'^(?P<uuid>' + uuidregex + ')/$', nrpDesigner, name = "nrpDesigner"),
					url(r'^(?P<uuid>' + uuidregex + ')/library/$', createLibrary, name="createLibrary"),
					url(r'^(?P<uuid>' + uuidregex + ')/library/process/$', processLibrary, name="processLibrary"),
					url(r'^(?P<uuid>' + uuidregex + ')/library/view/$', viewLibrary, name="viewLibrary"),
					url(r'^(?P<uuid>' + uuidregex + ')/library/download/$', downloadLibrary, name="downloadLibrary"),
					url(r'^nrpselection/(?P<uuid>' + uuidregex + ')/$', SpeciesListView, name="guiTool"),
					url(r'^selectableMonomers/$', get_available_monomers, name="selectableMonomers"),
					url(r'^makeConstruct/(?P<uuid>' + uuidregex + ')$', makeConstruct, name="makeConstruct"),
					url(r'^getConstruct/(?P<uuid>' + uuidregex + ')$', getConstruct, name="getConstruct"),
					url(r'^structure/', make_structure, name="Structure"),
					url(r'^submit/', submit_nrp, name="submitNRP"),
					url(r'^api/(?P<uuid>' + uuidregex + ')$', saveNrpMonomers, name = "saveNrpMonomers" ),
					url(r'^(?P<uuid>' + uuidregex + ')/domainSequence/', domainSequenceView, name="domainSequence")
					)
