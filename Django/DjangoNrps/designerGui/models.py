from django.db import models

import subprocess
import fcntl
import os
import select
from xml.dom.minidom import parseString
import json
import logging

from django.conf import settings
from django.db import connection
from databaseInput.models import Substrate, Domain
from gibson.models import Construct, ConstructFragment
from fragment.models import Gene

from celery.contrib.methods import task

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

    '''
Depiction of basic NRP model dependencies (ForeignKeys, ManyToManyFields):

 designerGui.NRP <----------------------- designerGui.DomainOrder  ---> databaseInput.Domain
   ||                                                       /\
   ||                                                       ||
   \/                                                       ||
  gibson.Construct <--- gibson.ConstructFragment ---> fragment.Gene


Some comments regarding deletion considerations/strategies:
*if domains get designed again(self.designDomains) -> delete (previous) DomainOrder [databaseInput.Domain objects obviously stay]
*if construct gets created again -> manually delete previous fragment.Gene objects (in case DomainOrder did not point yet
at appropriate fragment.Gene)-> gibson.ConstructFragment objects get deleted automatically
    '''

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

    def delete_dependencies(self):
        [x.gene.delete() for x in self.domainOrder.all() if x.gene is not None]
        self.designed = False
        self.domainOrder.all().delete()
        if self.construct is not None:
            [x.fragment.delete() for x in self.construct.cf.all() if x.fragment is not None]
            prev_construct = self.construct
            self.construct = None
            self.save()
            prev_construct.delete()
        pass

    def fullDelete(self):
        for gene in [x.gene for x in self.domainOrder.all()]:
            if gene is not None:
                gene.delete()
        if self.construct is not None:
            if len(self.construct.cf.all()) > 0 :
                [x.fragment.delete() for x in self.construct.cf.all()]
            self.construct.delete() #causes self to be deleted as well due to on cascade deletion in sql
        else:
            self.delete()

    # for non-celery use , e.g. via Shell, make sure that NRP has been designed before using method below
    def makeConstruct(self):
        if not self.designed:
            task = self.designDomains.delay()
            return task.id
        #import pdb;pdb.set_trace()
        # clean up previous stuff
        if self.construct is not None:
            if len(self.construct.cf.all()) > 0 :
                [x.fragment.delete() for x in self.construct.cf.all() if x.fragment is not None]
            prev_construct = self.construct
            self.construct = None
            self.save()
            prev_construct.delete()

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
                description = ' '.join(['DNA sequence of',
                    str(domain.domainType),
                    'domain of module',
                    str(domain.module),
                    'of gene',
                    str(domain.cds.geneName),
                    'of origin',
                    str(domain.cds.origin)]),
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
        return True

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

    @task()
    def designDomains(self):
        #see deletion strategy in multiline comment above
        if len(self.domainOrder.all())>0:
            [x.gene.delete() for x in self.domainOrder.all() if x.gene is not None]
            self.domainOrder.all().delete() #possibly redundant?
        self.designed = False
        lasterror = [None] # nonlocal only in python3
        logger = logging.getLogger('user_visible')
        xmlout = [""] # nonlocal only in python3
        def processErr(lines):
            lines = lines.splitlines()
            for line in lines:
                if line.startswith("WARNING: "):
                    logger.warning(line[9:])
                elif line.startswith("ERROR: "):
                    lasterror[0] = line[7:]
                    logger.error(lasterror[0])
                elif line.startswith("INFO: "):
                    logger.info(line[7:])
                else:
                    logger.error(line)
        def processOut(lines):
            xmlout[0] += lines
        # call NRPS Designer C++ program
        dbSettings = {
            '--mysql-host'     : settings.DATABASES[connection.alias]['HOST'],
            '--mysql-port'     : settings.DATABASES[connection.alias]['PORT'],
            '--mysql-user'     : settings.DATABASES[connection.alias]['USER'],
            '--mysql-password' : settings.DATABASES[connection.alias]['PASSWORD']
            }
        args = ["nrpsdesigner", "-m", self.getPeptideSequenceAsString()]
        for k,v in dbSettings.items():
            if v:
                args.extend([k, v])
        child = subprocess.Popen(args, 1, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

        fds = (child.stdout.fileno(), child.stderr.fileno())
        mask = fcntl.fcntl(fds[0], fcntl.F_GETFL)
        fcntl.fcntl(fds[0], fcntl.F_SETFL, mask | os.O_NONBLOCK)
        mask = fcntl.fcntl(fds[1], fcntl.F_GETFL)
        fcntl.fcntl(fds[1], fcntl.F_SETFL, mask | os.O_NONBLOCK)
        poll = select.poll()
        mask = select.POLLIN | select. POLLPRI | select.POLLERR | select.POLLHUP
        poll.register(fds[0], mask)
        poll.register(fds[1], mask)
        while child.poll() is None:
            fdsready = poll.poll()
            for fd in fdsready:
                if fd[0] == fds[1] and (fd[1] == select.POLLIN or fd[1] == select.POLLPRI):
                    readFrom(child.stderr, processErr)
                elif fd[0] == fds[0] and (fd[1] == select.POLLIN or fd[1] == select.POLLPRI):
                    readFrom(child.stdout, processOut)
        if child.returncode != 0:
            logger.error("NRPSDesigner returned errorcode %d" % child.returncode)
            raise Exception("NRPSDesigner returned errorcode %d: %s" % (child.returncode, lasterror[0]))

        # parse xml to extract domain list
        designerDom = parseString(xmlout[0])
        domainDomList = designerDom.getElementsByTagName('domain')
        domainIdList = [int(domainDom.getElementsByTagName('id')[0].firstChild.data) for domainDom in domainDomList]

        prevDomains = DomainOrder.objects.filter(nrp = self)
        prevDomains.delete()

        for count, domainId in enumerate(domainIdList):
            domain = Domain.objects.get(pk = domainId)
            domainOrder = DomainOrder.objects.create(nrp= self, domain=domain, order=count)
        self.designed = True
        self.save()

def readFrom(file, callback):
    try:
        lines = file.read()
        while len(lines) > 0:
            callback(lines)
            lines = file.read()
    except:
        pass

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


