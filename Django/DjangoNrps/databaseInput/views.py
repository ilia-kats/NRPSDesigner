from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.views.generic.base import TemplateView
from django.views.generic import CreateView
from django.core.urlresolvers import reverse_lazy
from django.template import loader, RequestContext

from django.views.generic.detail import DetailView
from django.views.generic.edit import CreateView
from django.contrib.auth import get_user_model

from django.forms.models import inlineformset_factory
import pdb

import json
import requests
import time
import pdb
from xml.dom.minidom import parseString
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from databaseInput.models import Origin, Cds, Domain
from databaseInput.forms import CdsFormSet, CdsForm, OriginForm, DomainForm

# Create your views here.

# class PfamView(TemplateView):
#     template_name = 'databaseInput/pfam.html'
#     #form_class = CdsForm(prefix = 'cds')
#     #success_url = reverse_lazy("pfam")
#     DomainFormSet = inlineformset_factory(Cds, Domain, extra=3)

#     def get_context_data(self, **kwargs):
#         context = super(PfamView, self).get_context_data(**kwargs)
#         context['form'] = CdsForm(prefix='cds')
#         context['originSet'] = DomainFormSet()
#         context['originForm'] = OriginForm(prefix='origin')
#         return context


#     def post(self,request,*args ,**kwargs):
#         if 'originSubmit' in request.POST:  #allow user to submit origin, if his is not already in db
#             originForm = OriginForm(request.POST, prefix = 'origin')
#             if originForm.is_valid():
#                 originForm.save()
#                 return HttpResponseRedirect(reverse_lazy("pfam"))
#             else:
#                 return HttpResponseRedirect(reverse_lazy("pfam"))
#         else:
#             cdsForm = CdsForm(request.POST, prefix= 'cds')
#             if cdsForm.is_valid():
#                 cds = cdsForm.save(commit=False)
#                 domainFormSet = DomainFormSet(request.POST, instance = cds)
#                 if domainFormSet.is_valid():
#                     cds.save()
#                     domainFormSet.save()
#                     return HttpResponseRedirect("http://www.google.com")
#                 else:
#                     return HttpResponseRedirect("http://www.yahoo.com")

#             return HttpResponseRedirect(reverse_lazy("pfam"))

class OriginCreateView(CreateView):
    template_name = 'databaseInput/addOrigin.html'
    form_class = OriginForm
    success_url=reverse_lazy('pfam')

def cdsInput(request):
    t = loader.get_template('databaseInput/cdsInput.html')
    
    c = RequestContext(request, {
        'form':CdsForm(prefix='cds'),
    })
    return HttpResponse(t.render(c))

def domainInput(request):
    if request.method == "POST":
        cdsId = int(request.POST['cdsId'])
        cds = Cds.objects.get(pk = cdsId)

        initialDict = cds.predictDomains()
        t = loader.get_template('databaseInput/domainInput.html')
        DomainFormSet = inlineformset_factory(Cds, Domain, form= DomainForm , extra=len(initialDict))
        c = RequestContext(request, {
        'originSet':DomainFormSet(initial = initialDict),
        })
        return HttpResponse(t.render(c))

class HomeTemplateView(TemplateView):
    template_name = 'home.html'

class ProfileTemplateView(TemplateView):
    template_name = 'databaseInput/profile.html'

class UserDetailView(DetailView):
    model = get_user_model()
    template_name = 'databaseInput/user_detail.html'
    slug_field = "username"

def sauceFunc(request):
    #first read FASTA file and translate sequence to protein!
    #dnaSeq  = SeqIO.read("bpsa.fasta", "fasta",IUPAC.unambiguous_dna).seq
    sequence = request.POST["sequence"].replace("\n","").replace("\t","")
    dnaSeq = Seq(sequence, IUPAC.unambiguous_dna)
    protSeq = dnaSeq.translate(to_stop=True)

    #send PFAM request 1
    pfamUrl = "http://pfam.sanger.ac.uk/search/sequence"
    pfamParams = {'seq':str(protSeq), 'output':'xml'}
    pfamRequest= requests.get(pfamUrl, params=pfamParams)

    #extract result link from XML by converting to DOM, extracting result url tag, removing tag elements from string
    pfamDom = parseString(pfamRequest.text)
    pfamJobId = pfamDom.getElementsByTagName('job')[0].attributes["job_id"].value
    pfamResultUrl = "http://pfam.sanger.ac.uk/search/sequence/resultset/" + pfamJobId
    pfamGraphicUrl = "http://pfam.sanger.ac.uk/search/sequence/graphic/" + pfamJobId

    #wait a bit
    time.sleep(1)

    # keep sending PFAM request 2 until it works
    while True:
        pfamResultRequest = requests.get(pfamResultUrl)
        if pfamResultRequest.status_code == 200:
            break
        elif pfamResultRequest.status_code == 202:
            time.sleep(1)
        else:
            break #THIS SHOULD actually throw exception!!!

    # read xml DOM, find stuff corresponding to domains
    #pfamResultDom = parseString(pfamResultRequest.text)
    #pfamResultMatches = pfamResultDom.getElementsByTagName('match')
    #pfamDomains = []
    #for match in pfamResultMatches:
    #   pfamDomains.append(match.attributes["id"].value)
    pfamResultRequest = requests.get(pfamGraphicUrl)
    import pdb; pdb.set_trace()
    return HttpResponse(pfamResultRequest.text[1:-1], mimetype="text/plain")

