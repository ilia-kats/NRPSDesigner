from django.conf.urls import patterns, url
from databaseInput.views import sauceFunc, cdsInput, domainInput, origin_add, origin_ajax_save, msa_domain_view


urlpatterns = patterns('',
    url(r'^submit/', sauceFunc, name="sauceFunc"),
    url(r'^$', cdsInput, name="pfam"),
    url(r'^addOrigin', origin_add, name="addOrigin"),
    url(r'^getDomainForm/', domainInput, name="getDomainForm"),
    url(r'^saveOrigin/', origin_ajax_save, name="saveOrigin"),
    url(r'^MSA/', msa_domain_view, name="msaDomainView")
)