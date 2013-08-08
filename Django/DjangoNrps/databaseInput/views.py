from django.shortcuts import render
from django.http import HttpResponse
from django.views.generic.base import TemplateView

# Create your views here.

class PfamView(TemplateView):
	template_name = 'databaseInput/pfam.html'
