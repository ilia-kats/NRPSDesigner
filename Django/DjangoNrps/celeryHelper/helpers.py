# completely failsafe script/functions!
# will be used to update celery task state,
# if there's any state to update or if 
# function was indeed called as celery task

from logging import Handler
#failsafe imports

try:
    from celery.result import AsyncResult
    from celery import current_task
except BaseException:
	pass

class TaskLogHandler(Handler):
    def emit(self, record):
        msg = self.format(record)
        print msg
        try:
            task = AsyncResult(current_task.request.id)
            print task
            try:
                log = task.result['log']
            except:
                log = []
            log.append(msg)
            current_task.update_state(state='log', meta={'log': log})
        except BaseException as e:
            pass
