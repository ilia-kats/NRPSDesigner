from django.shortcuts import render
from django.http import HttpResponse
from django.views.generic.base import TemplateView
from django.views.generic import CreateView
from django.core.urlresolvers import reverse_lazy

from django.views.generic.detail import DetailView
from django.contrib.auth import get_user_model
from databaseInput.models import Origin, Cds
from databaseInput.forms import CdsFormSet, DomainFormSet, CdsForm
# Create your views here.

class PfamView(CreateView):
	template_name = 'databaseInput/pfam.html'
	model = Cds
	success_url = reverse_lazy("pfam")

	def get_context_data(self, **kwargs):
		context = super(PfamView, self).get_context_data(**kwargs)
		context['originSet'] = DomainFormSet()
		return context

class HomeTemplateView(TemplateView):
	template_name = 'home.html'

class ProfileTemplateView(TemplateView):
	template_name = 'databaseInput/profile.html'

class UserDetailView(DetailView):
	model = get_user_model()
	template_name = 'databaseInput/user_detail.html'
	slug_field = "username"
