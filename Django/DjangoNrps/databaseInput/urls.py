from django.conf.urls import patterns, url
from databaseInput.views import sauceFunc, cdsInput, OriginCreateView, domainInput


urlpatterns = patterns('',
    url(r'^submit/', sauceFunc, name="sauceFunc"),
    url(r'^$', cdsInput, name="pfam"),
    url(r'^addOrigin', OriginCreateView.as_view(), name="addOrigin"),
    url(r'^getDomainForm/', domainInput, name="getDomainForm")
)