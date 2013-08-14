from django.conf.urls import patterns, include, url
from django.contrib import admin
from databaseInput.views import PfamView, UserDetailView, HomeTemplateView, ProfileTemplateView
from dajaxice.core import dajaxice_autodiscover, dajaxice_config
from django.contrib.staticfiles.urls import staticfiles_urlpatterns
from designerGui.views import SpeciesListView

dajaxice_autodiscover()
admin.autodiscover()

urlpatterns = patterns('',
    # Examples:
    # url(r'^$', 'DjangoNrps.views.home', name='home'),
    # url(r'^blog/', include('blog.urls')),
	url(r'^$', HomeTemplateView.as_view()),
	url(r'^admin/', include(admin.site.urls)),
    url(dajaxice_config.dajaxice_url, include('dajaxice.urls')),
    url(r'^pfam/', PfamView.as_view(), name="pfam"),
    url(r'^user/profile', ProfileTemplateView.as_view()),
    url(r'^user/', include('registration.backends.simple.urls')),
    url(r'^user/(?P<slug>\w+)/$', UserDetailView.as_view(), name="profile"),
    url(r'^tool/', SpeciesListView.as_view(), name="guiTool"),
)

urlpatterns += staticfiles_urlpatterns()
