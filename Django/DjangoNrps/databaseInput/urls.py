from django.conf.urls import patterns, url, include
from databaseInput.views import (cds_input, domain_prediction,
	origin_add, origin_ajax_save, msa_domain_view, do_msa,
	product_add, product_ajax_save, substrate_ajax_save,
    get_predicted_domain_formset,
	get_predicted_domain_formset_base, save_cds_domains,
    substrate_add, experiments)

# API patterns below used for entry of experimental evidence!
apipatterns = patterns('databaseInput.api',
    url(r'^getOrigins/$', 'get_origins', name="getOrigins"),
    url(r'^getDomainTypes/$', 'get_domain_types', name="getDomainTypes"),
    url(r'^getDomains/$', 'get_domains', name="getDomains"),
    url(r'^getDomainBioJsSequence/$' ,'get_domain_biojs_sequence', 
        name="getDomainBioJsSequence"),
)

urlpatterns = patterns('',
    url(r'^$', cds_input, name="pfam"),
    url(r'^addOrigin', origin_add, name="addOrigin"),
    url(r'^addProduct', product_add, name="addProduct"),
    url(r'^addSubstrate', substrate_add, name="addSubstrate"),
    url(r'^domainPrediction/', domain_prediction, name="startDomainPrediction"),
    url(r'^getDomainForm/(?P<task_id>\S+)', get_predicted_domain_formset, name="getPredictedDomainFormset"),
    url(r'^getDomainForm/', get_predicted_domain_formset_base, name="getPredictedDomainFormsetBase"),
    url(r'^saveCdsDomains/', save_cds_domains, name="saveCdsDomains"),
    url(r'^saveOrigin/', origin_ajax_save, name="saveOrigin"),
    url(r'^saveProduct/', product_ajax_save, name="saveProduct"),
    url(r'^saveSubstrate/', substrate_ajax_save, name="saveSubstrate"),
    url(r'^MSA/(?P<task_id>\S+)', msa_domain_view, name="msaDomainView"),
    url(r'^MSA/', do_msa, name="doMsa"),
    url(r'experiments/', experiments, name="experiments"),
    (r'^api/', include(apipatterns)),
)


