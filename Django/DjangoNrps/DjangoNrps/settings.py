"""
Django settings for DjangoNrps project.

For more information on this file, see
https://docs.djangoproject.com/en/dev/topics/settings/

For the full list of settings and their values, see
https://docs.djangoproject.com/en/dev/ref/settings/
"""

# Build paths inside the project like this: os.path.join(BASE_DIR, ...)
import os
BASE_DIR = os.path.dirname(os.path.dirname(__file__))

# django celery
import djcelery
djcelery.setup_loader()
BROKER_URL = 'django://'
CELERY_IMPORTS=("databaseInput.models", "designerGui.models")
CELERY_TRACK_STARTED = True
#CELERY_RESULT_BACKEND = "amqp"

# Quick-start development settings - unsuitable for production
# See https://docs.djangoproject.com/en/dev/howto/deployment/checklist/

# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = 'vu4ovjwz7gdljb&4l+!isd-dple+%tm337px02pg6720%g08oi'

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = True

TEMPLATE_DEBUG = True

ALLOWED_HOSTS = []

LOGIN_URL = "auth_login"

# Application definition

INSTALLED_APPS = (
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'databaseInput',
    'DjangoNrps', # for templatetags
    'registration',
    'designerGui',
    #'django_evolution',
    'annoying',
    'fragment',
    'gibson',
    'south',
    'djcelery',
    'celeryHelper',
    'kombu.transport.django',
)

MIDDLEWARE_CLASSES = (
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
)

ROOT_URLCONF = 'DjangoNrps.urls'

WSGI_APPLICATION = 'DjangoNrps.wsgi.application'


# Database
# https://docs.djangoproject.com/en/dev/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.mysql',
       #'NAME': os.path.join(BASE_DIR, 'nrps_designer'),
        'NAME': 'nrps_designer',
        'USER': 'root',
        'PASSWORD': '',
        'HOST': '127.0.0.1',
    }
}

# Internationalization
# https://docs.djangoproject.com/en/dev/topics/i18n/

LANGUAGE_CODE = 'en-us'

TIME_ZONE = 'UTC'

USE_I18N = True

USE_L10N = True

USE_TZ = True


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/dev/howto/static-files/
STATIC_URL = '/static/'

MEDIA_URL = '/media/'
MEDIA_ROOT = '/tmp/'

STATICFILES_DIRS = (
    os.path.join(BASE_DIR, "static"),
)

TEMPLATE_DIRS = (
   os.path.join(BASE_DIR, 'templates'),
)

TEMPLATE_LOADERS = (
    'django.template.loaders.filesystem.Loader',
    'django.template.loaders.app_directories.Loader',
    'django.template.loaders.eggs.Loader',
)

TEMPLATE_CONTEXT_PROCESSORS = (
    'django.contrib.auth.context_processors.auth',
    'django.core.context_processors.debug',
    'django.core.context_processors.i18n',
    'django.core.context_processors.media',
    'django.core.context_processors.static',
    'django.core.context_processors.request',
    'django.contrib.messages.context_processors.messages'
)

STATICFILES_FINDERS = (
    'django.contrib.staticfiles.finders.FileSystemFinder',
    'django.contrib.staticfiles.finders.AppDirectoriesFinder',
)

LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'verbose': {
            'format': '%(levelname)s %(asctime)s %(module)s %(process)d %(thread)d %(message)s'
        },
        'simple': {
            'format': '%(levelname)s %(message)s'
        },
    },
    'handlers': {
        'task_handler': {
            'level': 'INFO',
            'class': 'celeryHelper.helpers.TaskLogHandler',
        }
    },
    'loggers': {
        'user_visible': {
            'handlers': ['task_handler'],
            'propagate': True,
            'level': 'INFO',
        }
    }
}

LOGIN_REDIRECT_URL = 'userprofile'
ACCOUNT_ACTIVATION_DAYS = 7

EMAIL_BACKEND = 'django.core.mail.backends.smtp.EmailBackend'

# Host for sending e-mail.
EMAIL_HOST = 'localhost'

# Port for sending e-mail.
EMAIL_PORT = 25

# Optional SMTP authentication information for EMAIL_HOST.
EMAIL_HOST_USER = ''
EMAIL_HOST_PASSWORD = ''
EMAIL_USE_TLS = False

#caching
CACHES = {
    'default': {
        'BACKEND': 'django.core.cache.backends.db.DatabaseCache',
        'LOCATION': 'open_babel_structs',
    }
}

DEFAULT_FROM_EMAIL = 'kats@stud.uni-heidelberg.de'
CURATION_REQUEST_RECIPIENTS = ['nikos.ignatiadis01@gmail.com', 'k.herbst@stud.uni-heidelberg.de', 'nilskurzawa@yahoo.de', 'kats@stud.uni-heidelberg.de']
CURATION_GROUP = "curator"

UNAFOLD_WD = '/tmp/'
HYBRID_SS_MIN_PATH = 'hybrid-ss-min'
HYBRID_MIN_PATH = 'hybrid-min'
BOXPLOT_NG_PATH = 'boxplot_ng'

PREDEFINED_PART_TYPES = ['Expression Plasmid', 'Plasmid Backbone', 'Promoter', 'RBS', 'Terminator']
