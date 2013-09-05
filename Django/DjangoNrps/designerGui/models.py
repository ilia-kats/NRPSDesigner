from django.db import models

import subprocess
from xml.dom.minidom import parseString

from databaseInput.models import Substrate, Domain
from gibson.models import Construct, ConstructFragment
from fragment.models import Gene

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

	def __unicode__(self):
		return self.name

	def makeConstruct(self):
		if (not self.designed):
			self.designDomains()
		name = self.name + ' Gibson Construct'
		nrpConstruct = Construct.objects.create(owner = self.owner, 
			name = name,
			description = 'NRPS designer',
			shape = 'c')
		for count, domainId in enumerate(self.getDomainSequence()):
			domain = Domain.objects.get(pk=domainId)
			cdsSequence = domain.cds.dnaSequence
			domainStart = domain.definedStart - 1
			domainStop  = domain.definedLinkerStop
			domainSequence = cdsSequence[domainStart:domainStop]

			domainGene = Gene.objects.create(owner = self.owner,
				name = 'type:' + str(domain.domainType) + ' id:' + str(domainId),
				description = 'NRPS designer',
				sequence = domainSequence,
				origin = 'ND',
				viewable = 'H')

			domainConstructFragment = ConstructFragment.objects.create(
				construct = nrpConstruct,
				fragment = domainGene,
				order = count,
				direction = 'f'
				)

		self.construct = nrpConstruct
		self.save()

	def getPeptideSequence(self):
		monomers = self.monomers.order_by('substrateOrder').all()
		orderedMonomerIds = [int(monomer.pk) for monomer in monomers]
		return orderedMonomerIds

	def getPeptideSequenceAsString(self):
		peptideSequenceAsString = ','.join(["%s" % monomerId for monomerId in self.getPeptideSequence()])
		return peptideSequenceAsString

	def getPeptideSequenceForStructView(self):
		return '?' + '&'.join(["as=%s" % monomerId for monomerId in self.getPeptideSequence()])

	def getDomainSequence(self):
		domains = self.designerDomains.order_by('domainOrder').all()
		orderedDomainIds = [int(domain.pk) for domain in domains]
		return orderedDomainIds

	def designDomains(self):
		# call NRPS Designer C++ program
		rawXmlOutput = subprocess.check_output('nrpsdesigner '+self.getPeptideSequenceAsString(), shell=True)

		# parse xml to extract domain list
		designerDom = parseString(rawXmlOutput)
		domainDomList = designerDom.getElementsByTagName('domain')
		domainIdList = [int(domainDom.getElementsByTagName('id')[0].firstChild.data) for domainDom in domainDomList]

		prevDomains = DomainOrder.objects.filter(nrp = self)
		prevDomains.delete()

		for count, domainId in enumerate(domainIdList):
			domain = Domain.objects.get(pk = domainId)
			domainOrder = DomainOrder.objects.create(nrp= self, domain=domain, order=count)
		self.designed = True
		self.save()


class SubstrateOrder(models.Model):
	nrp = models.ForeignKey('NRP', related_name='substrateOrder')
	substrate = models.ForeignKey('databaseInput.Substrate', related_name='substrateOrder')
	order = models.PositiveIntegerField()

class DomainOrder(models.Model):
	nrp = models.ForeignKey('NRP', related_name = 'domainOrder')
	domain = models.ForeignKey('databaseInput.Domain', related_name = 'domainOrder')
	order = models.PositiveIntegerField()
	designerStart = models.IntegerField(blank=True, null=True)
	designerStop = models.IntegerField(blank=True, null=True)


