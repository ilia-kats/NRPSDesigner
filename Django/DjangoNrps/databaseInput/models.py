# This is an auto-generated Django model module.
# You'll have to do the following manually to clean this up:
#   * Rearrange models' order
#   * Make sure each model has one field with primary_key=True
#   * Remove `managed = False` lines for those models you wish to give write DB access
# Feel free to rename the models, but don't rename db_table values or field names.
#
# Also note: You'll have to insert the output of 'django-admin.py sqlcustom [appname]'
# into your database.
from __future__ import unicode_literals

from django.db import models
from django.contrib.contenttypes.models import ContentType
from django.contrib.contenttypes import generic



class Origin(models.Model):
    source = models.CharField(max_length=20) #usually taxon ID, can also be biobrick or plasmid identifier
    species = models.CharField(max_length = 100, blank=True, null= True)
    product = models.CharField(max_length=20)
    dnaSequence = models.TextField()
    description = models.TextField()
    #linkout = models.TextField()
    #norine = models.CharField(max_length=30, blank=True, null=True)
    linkout =  generic.GenericRelation('Linkout')
    parent = models.ForeignKey('self', blank=True, null=True, related_name='child')

    def __unicode__(self):
        return self.description

class Domain(models.Model):
    module = models.IntegerField()
    origin = models.ForeignKey('Origin')
    domainType = models.ForeignKey('Type')
    substrateSpecificity = models.ManyToManyField('Substrate', blank=True, null=True)
    chirality = models.CharField(max_length=1, choices= (('L','L'),('D','D')), blank=True, null=True)
    description = models.TextField()
    pfamLinkerStart = models.IntegerField()
    pfamLinkerStop = models.IntegerField()
    definedLinkerStart = models.IntegerField()
    definedLinkerStop = models.IntegerField()
    pfamStart = models.IntegerField()
    pfamStop = models.IntegerField()
    definedStart = models.IntegerField()
    definedStop = models.IntegerField()
    linkout =  generic.GenericRelation('Linkout')

    #linkout = models.TextField(blank=True, null=True)

    def __unicode__(self):
        return self.description

class Substrate(models.Model):
    name = models.CharField(max_length=30)
    chirality = models.CharField(max_length=1, choices= (('L','L'),('D','D')))
    structure = models.TextField()
    linkout =  generic.GenericRelation('Linkout')

    #linkout = models.TextField(blank=True, null=True)

    def __unicode__(self):
        return self.name

class Type(models.Model):
    name = models.CharField(max_length=4)
    description = models.TextField()

    def __unicode__(self):
        return self.name

class Linkout(models.Model):
    linkoutType = models.ForeignKey('LinkoutType')
    identifier = models.CharField(max_length=50)
    limit = models.Q(app_label = 'databaseInput', model = 'Substrate') | models.Q(app_label = 'databaseInput', model = 'Domain') | models.Q(app_label = 'databaseInput', model = 'Origin')
    content_type = models.ForeignKey(ContentType, limit_choices_to = limit)
    object_id = models.PositiveIntegerField()
    content_object = generic.GenericForeignKey('content_type', 'object_id')

    def __unicode__(self):
        return self.identifier

class LinkoutType(models.Model):
    shortcut = models.CharField(max_length=10)
    url = models.CharField(max_length=100)
    description = models.TextField(blank=True, null=True)

    def __unicode__(self):
        return self.shortcut
