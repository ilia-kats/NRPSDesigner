from django.db import models

import subprocess
from xml.dom.minidom import parseString
import json

from DjangoNrps.settings import DATABASES
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

	def fullDelete(self):
		for gene in [x.gene for x in self.domainOrder.all()]:
			if gene != None:
				ConstructFragment.objects.filter(fragment = gene).delete()
				gene.delete()
		self.construct.delete()
		self.delete()



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
			domainSequence = domain.get_sequence()

			domainGene = Gene.objects.create(owner = self.owner,
				name = 'type:' + str(domain.domainType),# + ' gene:' + str(domain.cds.geneName),
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

			domainOrder = DomainOrder.objects.get(nrp = self, domain=domain, order = count )
			domainOrder.gene = domainGene
			domainOrder.save()

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

	def generatePfamGraphicJson(self):
		domainList = self.getDomainSequence()
		domain_origins = []
		graphic_length = 350*len(domainList)
		regions = []
		i = 1
		for did in domainList:
			start = i*100
			end = start+200
			i += 3
			x_domain = Domain.objects.get(pk=did)
			domain_origins.append(x_domain.cds.origin)
			region_def = json.loads(x_domain.domainType.pfamGraphic)
			region_def.update({"start" : str(start), "end" : str(end)})
			regions.append(region_def)
		pfamJson = json.dumps({"length" : graphic_length, "regions": regions})
		return pfamJson


	def designDomains(self):
		# call NRPS Designer C++ program
		dbSettings = {
			'--mysql-host '     : DATABASES['default']['HOST'],
			'--mysql-port '     : DATABASES['default']['PORT'],
			'--mysql-user '     : DATABASES['default']['USER'],
			'--mysql-password  ': DATABASES['default']['PASSWORD']
			}
		dbCmds = ' '.join([k+v for k,v in dbSettings.items() if v!=''])
		rawXmlOutput = subprocess.check_output('nrpsdesigner ' +dbCmds +' -m '+self.getPeptideSequenceAsString(), shell=True)

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
	gene = models.ForeignKey('fragment.gene', null=True, blank=True)


