from django.conf.urls import patterns, url

from DjangoNrps import views

urlpatterns = patterns('',
    url(r'^$', views.index, name='index')
)