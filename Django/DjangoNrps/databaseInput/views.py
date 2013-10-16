from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseServerError, HttpResponseBadRequest
from django.views.generic.base import TemplateView
from django.views.generic import CreateView
from django.core.urlresolvers import reverse_lazy
from django.core.mail import send_mail
from django.template import loader, RequestContext

from django.views.generic.detail import DetailView
from django.views.generic.edit import CreateView
from django.contrib.auth import get_user_model

from django.forms.models import inlineformset_factory
from django.contrib.auth.decorators import login_required, user_passes_test
from django.conf import settings
from django.contrib.sites.models import Site
from django.contrib.sites.models import RequestSite

from celery.result import AsyncResult

import json
#import requests
import time

from xml.dom.minidom import parseString
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

from databaseInput.models import Origin, Cds, Domain
from databaseInput.forms import CdsFormSet, CdsForm, OriginForm, DomainForm, ProductForm

from gibson.jsonresponses import JsonResponse, ERROR

def do_msa(request):
    if request.method == "POST":
        cdsForm = CdsForm(request.POST, prefix="cds")
        #import pdb;pdb.set_trace()
        if cdsForm.is_valid():
            cds = cdsForm.save()
            # let's modify request.POST, so that the user still gets a
            # response, even if module number has not been entered yet
            post = request.POST.copy()
            post[u'domains-module'] = u'0'

            # prefix for domainForm extracted from DomainFormSet by JS
            domainForm = DomainForm(post,  prefix="domains")
            if domainForm.is_valid():

                initialDict = domainForm.cleaned_data

                del initialDict['substrateSpecificity']
                initialDict['cds'] = cds
                domain = Domain.objects.create(**initialDict)
                task = domain.align_same_type.delay()
                return JsonResponse({'taskId': task.id})
                t = loader.get_template('databaseInput/MSA_test.html')
                c = RequestContext(request,{
                    'jsonMSA': MSA
                })

                domain.delete(  )
                cds.delete()
                return HttpResponse(t.render(c))
        else:
            return HttpResponse("")      #think of how to best handle errors..
    else:
        return HttpResponse("")

def msa_domain_view(request, task_id):
    msa = AsyncResult(task_id).get()
    t = loader.get_template('databaseInput/MSA_test.html')
    c = RequestContext(request,{
        'jsonMSA': msa[1]
        })
    domain = Domain.objects.get(pk=msa[0])
    domain.cds.delete()
    domain.delete
    return HttpResponse(t.render(c))


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
#                     domainFormSet.save()pd
#                     return HttpResponseRedirect("http://www.google.com")
#                 else:
#                     return HttpResponseRedirect("http://www.yahoo.com")

#             return HttpResponseRedirect(reverse_lazy("pfam"))

@login_required
def product_add(request):
    t = loader.get_template("databaseInput/addProduct.html")
    c = RequestContext(request, {
        'form': ProductForm(prefix='product')
        })
    return HttpResponse(t.render(c))

@login_required
def origin_add(request):
    t = loader.get_template('databaseInput/addOrigin.html')
    c = RequestContext(request, {
        'form':OriginForm(prefix='origin'),
    })
    return HttpResponse(t.render(c))

@login_required
def product_ajax_save(request):
    if request.method == "POST":
        productForm = ProductForm(request.POST, prefix='product')
        if productForm.is_valid():
            product = productForm.save()
            product.user = request.user
            product.save()
            cdsForm = CdsForm(request.POST, prefix='cds')

            cdsForm.full_clean()
            initialDict = cdsForm.cleaned_data

            initialDict['product'] = product
            updatedCdsForm = CdsForm(initial=initialDict,prefix='cds')

            t = loader.get_template('databaseInput/cdsInputTab.html')
            c = RequestContext(request, {
            'form':updatedCdsForm,
                })
            return JsonResponse({'html': t.render(c)})

        else:
            t = loader.get_template('databaseInput/addProduct.html')
            c = RequestContext(request, {
            'form':productForm,
                })
            return JsonResponse({'html': t.render(c)}, ERROR)

@login_required
def origin_ajax_save(request):
    if request.method == "POST":
        originForm = OriginForm(request.POST, prefix='origin')
        if originForm.is_valid():
            origin = originForm.save()
            origin.user = request.user
            origin.save()
            cdsForm = CdsForm(request.POST, prefix='cds')

            cdsForm.full_clean()
            initialDict = cdsForm.cleaned_data

            initialDict['origin'] = origin
            updatedCdsForm = CdsForm(initial=initialDict,prefix='cds')

            t = loader.get_template('databaseInput/cdsInputTab.html')
            c = RequestContext(request, {
            'form':updatedCdsForm,
                })
            return JsonResponse({'html': t.render(c)})

        else:
            t = loader.get_template('databaseInput/addOrigin.html')
            c = RequestContext(request, {
            'form':originForm,
                })
            return JsonResponse({'html': t.render(c)}, ERROR)

@login_required
def cds_input(request):
    t = loader.get_template('databaseInput/cdsInput.html')
    cdsForm = CdsForm(prefix='cds')
    c = RequestContext(request, {
        'form':cdsForm,
        'isAjax':False
    })
    return HttpResponse(t.render(c))

@login_required
def domain_prediction(request):
    if request.method == "POST":
        cdsForm = CdsForm(request.POST, prefix='cds')
        if cdsForm.is_valid():
            cds = cdsForm.save(commit=False)
            #test_task = cds.predictDomains.delay()
            task = cds.predictDomains.delay()
            return JsonResponse({'taskId': task.id})
        else:
            t = loader.get_template('databaseInput/cdsInputTab.html')
            c = RequestContext(request, {
            'form':cdsForm,
            'isAjax': True,
                })
            return JsonResponse({'html': t.render(c)}, ERROR)

def get_predicted_domain_formset(request, task_id):
    if request.method == "POST":
        initialDict = AsyncResult(task_id).get()
        t = loader.get_template('databaseInput/domainInput.html')
        DomainFormSet = inlineformset_factory(Cds, Domain, form= DomainForm , extra=len(initialDict))
        c = RequestContext(request, {
        'domainFormSet':DomainFormSet(initial = initialDict),
            })
        #return HttpResponse(t.render(c))
        return JsonResponse({'html': t.render(c),})

# kind of hacky view in order to be able to call get_predicted_domain_formset
# with task_id from JS and be able to use reverse url lookup
def get_predicted_domain_formset_base(request):
    return HttpResponse("igem <3")

@login_required
def save_cds_domains(request):
    if request.method == "POST":
        DomainFormSet = inlineformset_factory(Cds, Domain, form= DomainForm)
        cdsForm = CdsForm(request.POST, prefix='cds')
        if cdsForm.is_valid():
            cds = cdsForm.save(commit=False)
            cds.user = request.user
            domainFormSet = DomainFormSet(request.POST, instance = cds)
            if domainFormSet.is_valid():
                cds.save()
                domains = domainFormSet.save(commit=False)
                for domain in domains:
                    domain.user = request.user
                    domain.save()
                domainFormSet.save_m2m()
                #if successful, render cds input form with clean form again!
                t = loader.get_template('databaseInput/cdsInputTab.html')
                cdsForm = CdsForm(prefix='cds')
                c = RequestContext(request, {
                        'form':cdsForm,
                        'isAjax':True,
                        'djangoSuccess': '<strong>Well done!</strong> Successfully added new coding sequence into database!'
                         })
                return JsonResponse({'html': t.render(c)})
            else:
                t = loader.get_template('databaseInput/domainInput.html')
                c = RequestContext(request, {
                     'domainFormSet':domainFormSet,
                    })
                return JsonResponse({'html': t.render(c)}, ERROR)
        else:
            # if cds form does not validate, return original page again..with appropriate error alert
            # this view is just in case and should not be returned under normal circumstances
            # as javascript checks for changes in input after prediction as well
            t = loader.get_template('databaseInput/cdsInputTab.html')
            c = RequestContext(request, {
                    'form':cdsForm,
                    'isAjax':True,
                    'djangoError': '<strong>Woops..</strong>It appears the original Cds got tampered with after domain prediction. Please fix your input and try again!'
                    })
            return JsonResponse({'html': t.render(c)})
    else:
        # should not be accessed via get request
        return HttpResponseBadRequest()

class HomeTemplateView(TemplateView):
    template_name = 'home.html'

class ProfileTemplateView(TemplateView):
    template_name = 'databaseInput/profile.html'

class UserDetailView(DetailView):
    model = get_user_model()
    template_name = 'databaseInput/user_detail.html'
    slug_field = "username"

@login_required
def change_password(request):
    if request.method == "POST" and "oldpw" in request.POST and "newpw" in request.POST:
        if not request.user.check_password(request.POST['oldpw']):
            return JsonResponse({"field": "oldpw"}, ERROR)
        else:
            request.user.set_password(request.POST['newpw'])
            request.user.save()
            return JsonResponse({})

@login_required
def request_curation_privs(request):
    if request.method == "POST" and "text" in request.POST:
        subject = loader.get_template('databaseInput/curation_email_subject.txt')
        text = loader.get_template('databaseInput/curation_email.txt')
        if Site._meta.installed:
            site = Site.objects.get_current()
        else:
            site = RequestSite(request)
        ctxt = RequestContext(request, {'request_text': request.POST['text'], 'site': site})
        subject = subject.render(ctxt)
        subject = ''.join(subject.splitlines())
        send_mail(subject, text.render(ctxt), request.user.email, settings.CURATION_REQUEST_RECIPIENTS)
    return HttpResponse()

#def sauceFunc(request):
    ##first read FASTA file and translate sequence to protein!
    ##dnaSeq  = SeqIO.read("bpsa.fasta", "fasta",IUPAC.unambiguous_dna).seq
    #sequence = request.POST["sequence"].replace("\n","").replace("\t","")
    #dnaSeq = Seq(sequence, IUPAC.unambiguous_dna)
    #protSeq = dnaSeq.translate(to_stop=True)

    ##send PFAM request 1
    #pfamUrl = "http://pfam.sanger.ac.uk/search/sequence"
    #pfamParams = {'seq':str(protSeq), 'output':'xml'}
    #pfamRequest= requests.get(pfamUrl, params=pfamParams)

    ##extract result link from XML by converting to DOM, extracting result url tag, removing tag elements from string
    #pfamDom = parseString(pfamRequest.text)
    #pfamJobId = pfamDom.getElementsByTagName('job')[0].attributes["job_id"].value
    #pfamResultUrl = "http://pfam.sanger.ac.uk/search/sequence/resultset/" + pfamJobId
    #pfamGraphicUrl = "http://pfam.sanger.ac.uk/search/sequence/graphic/" + pfamJobId

    ##wait a bit
    #time.sleep(1)

    ## keep sending PFAM request 2 until it works
    #while True:
        #pfamResultRequest = requests.get(pfamResultUrl)
        #if pfamResultRequest.status_code == 200:
            #break
        #elif pfamResultRequest.status_code == 202:
            #time.sleep(1)
        #else:
            #break #THIS SHOULD actually throw exception!!!

    ## read xml DOM, find stuff corresponding to domains
    ##pfamResultDom = parseString(pfamResultRequest.text)
    ##pfamResultMatches = pfamResultDom.getElementsByTagName('match')
    ##pfamDomains = []
    ##for match in pfamResultMatches:
    ##   pfamDomains.append(match.attributes["id"].value)
    #pfamResultRequest = requests.get(pfamGraphicUrl)
    #import pdb; pdb.set_trace()
    #return HttpResponse(pfamResultRequest.text[1:-1], mimetype="text/plain")
