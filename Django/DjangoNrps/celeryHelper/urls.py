from django.conf.urls import patterns, url
from celeryHelper.views import celery_task_log, celery_log_base


urlpatterns = patterns('',
	url(r'^log/(?P<task_id>\S+)', celery_task_log, name="celery_task_log"),
	url(r'^log/$', celery_log_base, name="celery_log_base")
	)





    
