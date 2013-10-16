from django.conf.urls import patterns, include, url
from django.contrib import admin
from databaseInput.views import  UserDetailView, HomeTemplateView, ProfileTemplateView, request_curation_privs, change_password
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from django.conf.urls.static import static
from django.contrib.auth import login
from django.conf import settings
from registration import signals

admin.autodiscover()

urlpatterns = patterns('',
    # Examples:
    # url(r'^$', 'DjangoNrps.views.home', name='home'),
    # url(r'^blog/', include('blog.urls')),
	url(r'^$', HomeTemplateView.as_view(), name='home'),
	url(r'^admin/', include(admin.site.urls)),
	url(r'^user/requestCurationPrivs', request_curation_privs, name="submitCurationRequest"),
	url(r'^user/changePassword', change_password, name="changePassword"),
    url(r'^user/profile', ProfileTemplateView.as_view(), name="userprofile"),
    url(r'^user/', include('registration.backends.default.urls')),
    url(r'^users/(?P<slug>\w+)/$', UserDetailView.as_view(), name="profile"),
    url(r'^tool/', include('designerGui.urls')),
    url(r'^gibthon/', include('gibson.urls')),
    url(r'^fragment/', include('fragment.urls')),
    url(r'^input/', include('databaseInput.urls')),
    url(r'^celery/',include('celeryHelper.urls')),
)

urlpatterns += staticfiles_urlpatterns()

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)

def login_on_activation(user, request, **kwargs):
    user.backend='django.contrib.auth.backends.ModelBackend'
    login(request, user)

signals.user_activated.connect(login_on_activation)
