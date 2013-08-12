from django.shortcuts import render
from django.http import HttpResponse
from django.views.generic.base import TemplateView
from django.views.generic import CreateView

from django.views.generic.detail import DetailView
from django.contrib.auth import get_user_model
from databaseInput.models import Origin
# Create your views here.

class PfamView(CreateView):
	template_name = 'databaseInput/pfam.html'
	model = Origin

class HomeTemplateView(TemplateView):
	template_name = 'home.html'

class ProfileTemplateView(TemplateView):
	template_name = 'databaseInput/profile.html'

class UserDetailView(DetailView):
	model = get_user_model()
	template_name = 'databaseInput/user_detail.html'
	slug_field = "username"
