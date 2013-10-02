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
   ||
   ||
   \/
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
    indigoidineTagged = models.BooleanField(default=False)
    designerDomains = models.ManyToManyField('databaseInput.Domain', through = 'DomainOrder', blank=True, related_name = 'includedIn')
    construct = models.ForeignKey('gibson.Construct', null=True, blank=True)

    def __unicode__(self):
        return self.name

    def delete_dependencies(self, all=True):
        self.designed = False
        self.domainOrder.all().delete()
        if all:
            SubstrateOrder.objects.filter(nrp = self).delete()
        if self.construct is not None:
            for x in self.construct.cf.all():
                if x.fragment is not None and x.fragment.origin == 'ND':
                    x.fragment.delete()
                x.delete()
            if self.construct.primer is not None:
                [x.del_all() for x in self.construct.primer.all() if x is not None]
            self.construct.reset()
        self.save()

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
        # helper function
        # check whether 2 domains are placed next to each other in the same coding sequence
        def tupleIsConnected((x,y)):
            if x.cds != y.cds:
                return False
            cdsDomainList = x.cds.get_ordered_domain_list()
            if cdsDomainList.index(y) == cdsDomainList.index(x) + 1:
                return True
            else:
                return False

        # actual function
        if not self.designed or len(self.domainOrder.all()) == 0:
            task = self.designDomains.delay()
            return task.id

        domains = list(self.getDomainModelSequence())
        connectedDomainList = list()

        connectedDomains = [domains[0]]
        for domainTuple in zip(domains,domains[1:]):
            if not tupleIsConnected(domainTuple):
                connectedDomainList.append(connectedDomains)
                connectedDomains = list()
            connectedDomains.append(domainTuple[1])
        connectedDomainList.append(connectedDomains)

        # start creating construct
        name = self.name + ' Gibson Construct'
        if self.construct is None:
            self.construct = Construct.objects.create(owner = self.owner,
                name = name,
                description = 'NRPS designer',
                shape = 'c')

        # each list of connectedDomains corresponds to 1 fragment.gene and hence to 1 construct fragment
        for count, connectedDomains in enumerate(connectedDomainList):

            domain1 = connectedDomains[0]
            domain2 = connectedDomains[-1]
            domainSequence = domain1.cds.get_sequence(domain1,domain2)

            domainGene = Gene.objects.create(owner = self.owner,
                name = ','.join([x.domainType.name for x in connectedDomains]),
                description = ' '.join(['DNA sequence of',
                        ','.join(map(lambda x: ' '.join([str(x.domainType),
                                            'domain of module',
                                            str(x.module)]),
                                        connectedDomains)),
                    'of gene',
                    str(domain1.cds.geneName),
                    'of origin',
                    str(domain1.cds.origin)]),
                sequence = domainSequence,
                origin = 'ND',
                viewable = 'H')

            domainConstructFragment = ConstructFragment.objects.create(
                construct = self.construct,
                fragment = domainGene,
                order = count,
                direction = 'f'
                )

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

    #returns actual domain objects
    def getDomainModelSequence(self):
        domains = self.designerDomains.order_by('domainOrder').all()
        return domains

    #returns IDs
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
    def designDomains(self, curatedonly=True):
        #see deletion strategy in multiline comment above
        self.delete_dependencies(False)
        lasterror = [None] # nonlocal only in python3
        logger = logging.getLogger('user_visible')
        xmlout = [""] # nonlocal only in python3
        def processErr(lines):
            cat = 0
            def log():
                if lasterror[0]is not None:
                    return
                if cat == 0:
                    logger.error(lasterror[0])
                elif cat == 1:
                    logger.warning(lasterror[0])
                elif cat == 2:
                    logger.info(lasterror[0])
            lines = lines.splitlines()
            for line in lines:
                if line.startswith("WARNING: "):
                    log()
                    cat = 1
                    lasterror[0] = line[9:]
                elif line.startswith("ERROR: "):
                    log()
                    cat = 0
                    lasterror[0] = line[7:]
                elif line.startswith("INFO: "):
                    log()
                    cat = 2
                    lasterror[0] = line[6:]
                elif line.startswith(" "):
                    lasterror[0] += line.strip()
                else:
                    log()
                    cat = 0
                    lasterror[0] = line
            log()
        def processOut(lines):
            xmlout[0] += lines
        # call NRPS Designer C++ program
        dbSettings = {
            '--mysql-host'     : settings.DATABASES[connection.alias]['HOST'],
            '--mysql-port'     : settings.DATABASES[connection.alias]['PORT'],
            '--mysql-user'     : settings.DATABASES[connection.alias]['USER'],
            '--mysql-password' : settings.DATABASES[connection.alias]['PASSWORD']
            }
        args = ["nrpsdesigner", "-m", self.getPeptideSequenceAsString(), "-s", "-"]
        if curatedonly:
            args.extend(["--curated-only", "--curation-group", settings.CURATION_GROUP])
        if self.indigoidineTagged:
            args.append("-t")
        for k,v in dbSettings.items():
            if v:
                args.extend([k, v])
        child = subprocess.Popen(args, 0, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

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
        if child.returncode == -11:
            lasterror[0] = "Segmentation fault"
        if child.returncode != 0:
            logger.error("NRPSDesigner returned errorcode %d" % child.returncode)
            raise Exception("NRPSDesigner returned errorcode %d: %s" % (child.returncode, lasterror[0]))

        # parse xml to extract domain list
        # can not use libsbol here, as it only supports reading from file
        designerDom = parseString(xmlout[0])
        domainDomList = designerDom.getElementsByTagName('s:DnaComponent')
        domainIdList = []
        for comp in domainDomList:
            cid = comp.attributes.getNamedItem("rdf:about").value[1:]
            if cid.isdigit():
                for didn in comp.getElementsByTagName('s:displayId'):
                    did = didn.firstChild.data
                    if did[0] == '_':
                        domainIdList.append(int(did[1:]))
        prevDomains = DomainOrder.objects.filter(nrp = self)
        prevDomains.delete()

        for count, domainId in enumerate(domainIdList):
            domain = Domain.objects.get(pk = domainId)
            domainOrder = DomainOrder.objects.create(nrp= self, domain=domain, order=count)
        self.designed = True
        self.save()
        self.makeConstruct()

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
