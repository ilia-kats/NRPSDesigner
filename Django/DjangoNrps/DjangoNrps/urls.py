from django.conf.urls import patterns, include, url
from django.contrib import admin
from databaseInput.views import  UserDetailView, HomeTemplateView, ProfileTemplateView
from django.contrib.staticfiles.urls import staticfiles_urlpatterns

admin.autodiscover()

urlpatterns = patterns('',
    # Examples:
    # url(r'^$', 'DjangoNrps.views.home', name='home'),
    # url(r'^blog/', include('blog.urls')),
	url(r'^$', HomeTemplateView.as_view()),
	url(r'^admin/', include(admin.site.urls)),
    url(r'^user/profile', ProfileTemplateView.as_view()),
    url(r'^user/', include('registration.backends.simple.urls')),
    url(r'^user/(?P<slug>\w+)/$', UserDetailView.as_view(), name="profile"),
    url(r'^tool/', include('designerGui.urls')),
    url(r'^gibthon/', include('gibson.urls')),
    url(r'^fragment/', include('fragment.urls')),
    url(r'^input/', include('databaseInput.urls')),
)

urlpatterns += staticfiles_urlpatterns()
