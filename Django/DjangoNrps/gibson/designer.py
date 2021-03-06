from django.contrib.auth.decorators import login_required
from django.template import Context, loader, RequestContext
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound, Http404
from views import get_construct
from django.conf import settings

from models import *
from forms import *

def designer(request,cid):
    con = get_construct(request.user, cid)
    if con:
        t = loader.get_template('gibson/designer.html')
        c = RequestContext(request,{
            #'title':'Construct Designer',
            #'id': cid,
            #'construct':con,
        })
        return HttpResponse(t.render(c))
    else:
        raise Http404()

def design_tab(request, cid):
    con = get_construct(request.user, cid)
    if con:
        t = loader.get_template('gibson/designtab.html')
        c = RequestContext(request,{
            #'id': cid,
            'construct':con,
            'part_types': settings.PREDEFINED_PART_TYPES
        })
        return HttpResponse(t.render(c))
    else:
        return HttpResponseNotFound()

def construct_settings(request, cid):
    con = get_construct(request.user, cid)
    if con:
        if request.method == 'POST':
            form = SettingsForm(request.POST, instance=con.settings)
            if form.is_valid():
                form.save()
                return HttpResponse()
        t = loader.get_template('gibson/settings.html')
        s = con.settings
        c = RequestContext(request, {
            'cid': cid,
            'settings':s,
        })
        return HttpResponse(t.render(c))
    else:
        return HttpResponseNotFound()
