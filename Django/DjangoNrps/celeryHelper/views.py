import json

from gibson.jsonresponses import JsonResponse, ERROR

from celery.result import AsyncResult
def celery_task_log(request, task_id):
    task = AsyncResult(task_id)
    out = dict()
    if task.ready() and task.successful():
        out = {'status':'SUCCESS'}
    elif task.ready() and task.failed():
        out = {'status': 'FAILED', 'output': str(task.result)}
    elif task.status == "log":
        if 'log' in task.result:
            task_log = task.result['log']
        else:
            task_log = []
        out = {'status':'log', 'output':task_log}
    else:
        out = {'status': task.status}
    out['taskId'] = task_id;
    return JsonResponse(out)

def celery_log_base(request):
    return HttpResponse("UGA UGA said the mRNA to the ribosome making it cry")