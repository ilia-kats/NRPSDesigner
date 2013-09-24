import json

from django.http import HttpResponse

from celery.result import AsyncResult
def celery_task_log(request, task_id):
	task = AsyncResult(task_id)
	if task.ready():
		json_log = json.dumps({'status':'ready'})
	elif task.status == "log":
		task_log = task.result['log']
		json_log = json.dumps({'status':'log', 'output':task_log})
	else:
		json_log = json.dumps({'status': task.status})
	return HttpResponse(json_log,  content_type="application/json")

def celery_log_base(request):
	return HttpResponse("UGA UGA said the mRNA to the ribosome making it cry")