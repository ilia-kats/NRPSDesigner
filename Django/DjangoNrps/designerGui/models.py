from django.db import models
from databaseInput.models import Substrate, Domain
from gibson.models import Construct

class Species(models.Model):
    species = models.CharField(max_length=100)
    taxon_id = models.CharField(max_length=20)
  
    def __unicode__(self):
        return str(species)

    class Meta:
        verbose_name = "Species"
        verbose_name_plural = "Species"
        verbose_name = "Modifications"
        verbose_name_plural = "Modifications"

class NRP(models.Model):
	owner = models.ForeignKey('auth.User', null=True)
	name = models.CharField(max_length=80)
	description = models.CharField(max_length=2000, null= True, blank=True)
	monomers = models.ManyToManyField('databaseInput.Substrate', through='SubstrateOrder', blank=True, related_name='includedIn')
	created = models.DateTimeField(auto_now_add=True)
	modified = models.DateTimeField(auto_now=True)
	designed = models.BooleanField(default=False)
	designerDomains = models.ManyToManyField('databaseInput.Domain', through = 'DomainOrder', blank=True, related_name = 'includedIn')
	construct = models.ForeignKey('gibson.Construct', null=True, blank=True)

class SubstrateOrder(models.Model):
	nrp = models.ForeignKey('NRP', related_name='substrateOrder')
	substrate = models.ForeignKey('databaseInput.Substrate', related_name='substrateOrder')
	order = models.PositiveIntegerField()

class DomainOrder(models.Model):
	nrp = models.ForeignKey('NRP', related_name = 'domainOrder')
	domain = models.ForeignKey('databaseInput.Domain', related_name = 'domainOrder')
	order = models.PositiveIntegerField()
	designerStart = models.IntegerField()
	designerStop = models.IntegerField()


