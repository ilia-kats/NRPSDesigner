from django.conf.urls import patterns, url
from databaseInput.views import sauceFunc, cdsInput, OriginCreateView


urlpatterns = patterns('',
    url(r'^submit/', sauceFunc, name="sauceFunc"),
    url(r'^$', cdsInput, name="pfam"),
    url(r'^addOrigin', OriginCreateView.as_view(), name="addOrigin")
)