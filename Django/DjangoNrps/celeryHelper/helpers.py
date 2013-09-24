# completely failsafe script/functions!
# will be used to update celery task state,
# if there's any state to update or if 
# function was indeed called as celery task

import logging
#failsafe imports

try:
    from celery.result import AsyncResult
    from celery import current_task
except BaseException:
	pass
# failsafe function which will be called at different points of execution

def update_celery_task_state_log(msg):
	try:
		task = AsyncResult(current_task.request.id)
		try:
			old_msg = task.result['log']
		except:
			old_msg = ""
		new_msg = old_msg  + msg + "<br/>"
		current_task.update_state(state='log',
			meta={'log':new_msg})
	except BaseException:
		pass
	logging.info(msg)