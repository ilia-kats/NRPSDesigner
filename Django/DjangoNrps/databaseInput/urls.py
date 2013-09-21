from django.conf.urls import patterns, url
from databaseInput.views import (sauceFunc, cds_input, domain_prediction, 
	origin_add, origin_ajax_save, msa_domain_view, 
	product_add, product_ajax_save)


urlpatterns = patterns('',
    url(r'^submit/', sauceFunc, name="sauceFunc"),
    url(r'^$', cds_input, name="pfam"),
    url(r'^addOrigin', origin_add, name="addOrigin"),
    url(r'^addProduct', product_add, name="addProduct"),
    url(r'^getDomainForm/', domain_prediction, name="getDomainForm"),
    url(r'^saveOrigin/', origin_ajax_save, name="saveOrigin"),
    url(r'^saveProduct/', product_ajax_save, name="saveProduct"),
    url(r'^MSA/', msa_domain_view, name="msaDomainView")
)