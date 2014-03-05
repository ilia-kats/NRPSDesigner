from django.db import models

import subprocess
import fcntl
import os
import select
from xml.dom.minidom import parseString, getDOMImplementation
import json
import logging
from uuid import uuid4

from django.conf import settings
from django.db import connection
from databaseInput.models import Substrate, Domain
from gibson.models import Construct, ConstructFragment, fragment_feature, Settings, PCRSettings
from fragment.models import Gene, DomainGene, Feature, Qualifier

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

def make_uuid():
    return str(uuid4())

class NRP(models.Model):
    owner = models.ForeignKey('auth.User', null=True)
    uuid = models.CharField(max_length=36, db_index=True, default=make_uuid)
    name = models.CharField(max_length=80)
    description = models.CharField(max_length=2000, null= True, blank=True)
    monomers = models.ManyToManyField('databaseInput.Substrate', through='SubstrateOrder', blank=True, related_name='includedIn')
    created = models.DateTimeField(auto_now_add=True)
    modified = models.DateTimeField(auto_now=True)
    designed = models.BooleanField(default=False)
    indigoidineTagged = models.BooleanField(default=False)
    designerDomains = models.ManyToManyField('databaseInput.Domain', through = 'DomainOrder', blank=True, related_name = 'includedIn')
    construct = models.OneToOneField('gibson.Construct', null=True, blank=True, related_name='nrp')
    parent = models.ForeignKey('self', blank=True, null=True, related_name='child')
    boundary_parent = models.ForeignKey('self', blank=True, null=True, related_name='boundary_child')

    class Meta:
        verbose_name = "Nonribosomal peptide"
        verbose_name_plural = "Nonribosomal peptides"

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
                    Gene.remove(x.fragment.owner, x.fragment.pk)
                x.delete()
            if self.construct.primer is not None:
                [x.del_all() for x in self.construct.primer.all() if x is not None]
            self.construct.reset()
        self.parent = None
        self.save()

    def fullDelete(self):
        if self.construct is not None:
            self.delete_dependencies()
            self.construct.delete() #causes self to be deleted as well due to on cascade deletion in sql
        else:
            self.delete()

    # for non-celery use , e.g. via Shell, make sure that NRP has been designed before using method below
    def makeConstruct(self):
        # helper function
        # check whether 2 domains are placed next to each other in the same coding sequence
        def tupleIsConnected((x,y)):
            if x.domain.cds != y.domain.cds:
                return False
            cdsDomainList = x.domain.cds.get_ordered_domain_list()
            adjacent = (cdsDomainList.index(y.domain) == cdsDomainList.index(x.domain) + 1)
            if y.left_boundary is None:
                ystart = y.domain.get_start(x.domain, with_linker=True)
            else:
                ystart = y.left_boundary
            if x.right_boundary is None:
                xstop = x.domain.get_stop(y.domain, with_linker=True)
            else:
                xstop = x.right_boundary
            continuous_seq = (ystart == xstop)
            if adjacent and (y.left_boundary is None and x.right_boundary is None or continuous_seq):
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
            if self.description:
                description = ':\n' + self.description
            else:
                description = ''
            self.construct = Construct.objects.create(owner = self.owner,
                name = name,
                description = 'NRPS designer' + description,
                shape = 'c')

        lastDomain = None
        nextDomain = None

        # each list of connectedDomains corresponds to 1 fragment.gene and hence to 1 construct fragment
        for count, connectedDomains in enumerate(connectedDomainList):
            domain1 = connectedDomains[0]
            domain2 = connectedDomains[-1]

            if count < len(connectedDomainList) - 1:
                nextDomain = connectedDomainList[count + 1][0].domain.pk
            else:
                nextDomain = None

            if domain1.left_boundary is None:
                start_pos = domain1.domain.get_start(lastDomain)
                seqstart = lastDomain
            else:
                start_pos = seqstart = domain1.left_boundary
            if domain2.right_boundary is None:
                seqstop = domain2.domain
            else:
                seqstop = domain2.right_boundary
            domainSequence = domain1.domain.cds.get_sequence(domain1.domain,domain2.domain, seqstart, seqstop)

            domainGene = DomainGene.objects.create(owner = self.owner,
                name = ','.join([x.domain.domainType.name for x in connectedDomains]),
                description = ' '.join(['DNA sequence of',
                        ','.join(map(lambda x: ' '.join([str(x.domain.domainType),
                                            'domain of module',
                                            str(x.domain.module)]),
                                        connectedDomains)),
                    'of gene',
                    str(domain1.domain.cds.geneName),
                    'of origin',
                    str(domain1.domain.cds.origin)]),
                sequence = domainSequence,
                origin = 'ND',
                viewable = 'H')
            for (i,domain) in enumerate(connectedDomains):
                domainGene.domains.add(domain.domain)
                if domain.left_boundary is None:
                    if i == 0:
                        lastd = lastDomain
                    else:
                        lastd = None
                    start = domain.domain.get_start(lastd, with_linker=False)
                else:
                    start = domain.left_boundary
                if domain.right_boundary is None:
                    if i == len(connectedDomains):
                        nextd = nextDomain
                    else:
                        nextd = None
                    stop = domain.domain.get_stop(nextd, with_linker=False)
                else:
                    stop = domain.right_boundary
                f = Feature(type="domain", start=start - start_pos, end=stop - start_pos, direction='f', gene=domainGene)
                f.save()
                Qualifier(name="type", data=domain.domain.domainType.name, feature=f).save()
                Qualifier(name="module", data=domain.domain.module, feature=f).save()
                if domain.domain.substrateSpecificity.count() > 0:
                    Qualifier(name="substrate", data=", ".join([s.name for s in domain.domain.substrateSpecificity.all()]), feature=f).save()
                curated = "No"
                if len(domain.domain.user.groups.filter(name=settings.CURATION_GROUP)) > 0:
                    curated = "Yes"
                Qualifier(name="curated", data=curated, feature=f).save()
            f = Feature(type=fragment_feature, start=0, end=len(domainSequence), direction='f', gene=domainGene)
            f.save()
            Qualifier(name="gene", data=domain1.domain.cds.geneName, feature=f).save()
            Qualifier(name="species", data=domain1.domain.cds.origin.species, feature=f).save()
            if domain1.domain.cds.origin.sourceType == 'Species':
                sname = 'taxID'
            elif domain1.domain.cds.origin.sourceType == 'Biobrick':
                sname = 'BioBrick'
            else:
                sname = domain1.domain.cds.origin.sourceType
            Qualifier(name=sname, data=domain1.domain.cds.origin.source, feature=f).save()

            domainConstructFragment = ConstructFragment.objects.create(
                construct = self.construct,
                fragment = domainGene,
                order = count,
                direction = 'f'
                )

            lastDomain = domain2.domain

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
        if self.pk:
            return '?' + '&'.join(["as=%s" % monomerId for monomerId in self.getPeptideSequence()])
        else:
            return '?'

    #returns actual domain objects
    def getDomainModelSequence(self):
        domains = DomainOrder.objects.filter(nrp=self).order_by('order')
        return domains

    #returns IDs
    def getDomainSequence(self):
        domains = DomainOrder.objects.filter(nrp=self).order_by('order')
        orderedDomainIds = [int(domain.domain.pk) for domain in domains]
        return orderedDomainIds

    def generatePfamGraphicJson(self):
        domainList = self.getDomainSequence()
        graphic_length = 350*len(domainList)
        regions = []
        i = 1
        for did in domainList:
            start = i*100
            end = start+200
            i += 3
            x_domain = Domain.objects.get(pk=did)
            substrates = []
            for substrate in x_domain.substrateSpecificity.all():
                substrates.append({'name': substrate.name, 'chirality': substrate.chirality})
            region_def = json.loads(x_domain.domainType.pfamGraphic)
            region_def.update({"start" : str(start), "end" : str(end), "metadata": {'type': x_domain.domainType.name, 'substrates': substrates, 'chirality': x_domain.chirality, 'description': x_domain.description, 'curated': len(x_domain.user.groups.filter(name=settings.CURATION_GROUP)) > 0, 'gene': x_domain.cds.geneName, 'species': x_domain.cds.origin.species, 'source': x_domain.cds.origin.source, 'sourceType': x_domain.cds.origin.sourceType}})
            regions.append(region_def)
        pfamJson = json.dumps({"length" : graphic_length, "regions": regions})
        return pfamJson

    def adjustConstruct(self, target):
        def porder(order, construct):
            porder = order - 1
            if porder < 0:
                porder = construct.cf.count() - 1
            return porder
        def norder(order, construct):
            norder = order + 1
            if norder >= construct.cf.count():
                norder = 0
            return norder
        if not self.construct.processed:
            for j in self.construct.process():
                pass
        for (md, (sobj, nobj)) in {Settings: (self.construct.settings, target.construct.settings), PCRSettings: (self.construct.pcrsettings, target.construct.pcrsettings)}.iteritems():
            for field in md._meta.fields:
                if not isinstance(field, models.ManyToManyField) and not isinstance(field, models.OneToOneField) and not isinstance(field, models.ForeignKey) and not isinstance(field, models.AutoField):
                    setattr(nobj, field.name, getattr(sobj, field.name))
            nobj.save()
        ncfs = {}
        for j,cf in enumerate(self.construct.cf.order_by('order')):
            ncf = None
            if cf.fragment.origin != 'ND':
                ncf = target.construct.add_fragment(cf.fragment, j, cf.direction)
            else:
                try:
                    isIdentical = True
                    if cf.fragment.domaingene.domains.count() == target.construct.cf.get(order=j).fragment.domaingene.domains.count():
                        ourdomains = cf.fragment.domaingene.domains.all()
                        otherdomains = target.construct.cf.get(order=j).fragment.domaingene.domains.all()
                        for i in xrange(len(ourdomains)):
                            if ourdomains[i].pk != otherdomains[i].pk:
                                isIdentical = False
                                break
                    else:
                        isIdentical = False
                    if isIdentical:
                        ncf = target.construct.cf.get(order=j)
                except DomainGene.DoesNotExist:
                    pass
            if ncf is not None:
                ncf.concentration = cf.concentration
                ncf.save()
                ncfs[cf] = ncf
        for j in target.construct.process():
            pass

        for (cf, ncf) in ncfs.items():
            ptop = ncf.primer_top()
            pbottom = ncf.primer_bottom()
            ptop.stick.length = cf.primer_top().stick.length
            pbottom.stick.length = cf.primer_bottom().stick.length
            pf = self.construct.cf.get(order=porder(cf.order, self.construct))
            nf = self.construct.cf.get(order=norder(cf.order, self.construct))
            if pf in ncfs and ncfs[pf].order == porder(ncf.order, ncf.construct):
                pbottom.flap.length = cf.primer_bottom().flap.length
            if nf in ncfs and ncfs[nf].order == norder(ncf.order, ncf.construct):
                ptop.flap.length = cf.primer_top().flap.length
            ptop.stick.save()
            pbottom.stick.save()
            ptop.flap.save()
            pbottom.flap.save()
            ncf.construct.reprocess_primer(ptop)
            ncf.construct.reprocess_primer(pbottom)

    @task()
    def designDomains(self, curatedonly=True):
        #see deletion strategy in multiline comment above
        self.delete_dependencies(False)
        args = ["-m", self.getPeptideSequenceAsString(), "-s", "-"]
        if curatedonly:
            args.extend(["--curated-only", "--curation-group", settings.CURATION_GROUP])
        if self.indigoidineTagged:
            args.append("-t")
        xmlout = self._runNrpsDesigner(args)
        # parse xml to extract domain list
        # can not use libsbol here, as it only supports reading from file
        designerDom = parseString(xmlout)

        domainIdList = self._parseNrps(designerDom)
        self._design(domainIdList)

    @task()
    def makeLibrary(self, monomers, curatedonly=True):
        self.child.all().delete()
        impl = getDOMImplementation()
        xml = impl.createDocument(None, 'nrp', None)
        root = xml.documentElement
        monomersnode = xml.createElement('monomers')
        varpos = -1
        i = 0
        for monomer in monomers:
            monomernode = xml.createElement('monomer')
            if len(monomer) > 1:
                varpos = i
            i += 1
            for aa in monomer:
                idnode = xml.createElement("id")
                aanode = xml.createTextNode(str(aa))
                idnode.appendChild(aanode)
                monomernode.appendChild(idnode)
            monomersnode.appendChild(monomernode)
        root.appendChild(monomersnode)
        nrpsnode = xml.createElement('nrps')
        domains = self.getDomainSequence()
        for domain in domains:
            domainnode = xml.createElement('domain')
            idnode = xml.createElement('id')
            idtext = xml.createTextNode(str(domain))
            idnode.appendChild(idtext)
            domainnode.appendChild(idnode)
            nrpsnode.appendChild(domainnode)
        root.appendChild(nrpsnode)
        xmlstr = xml.toxml("utf-8")
        args = ['-n', '-', '-s', '-']
        if curatedonly:
            args.extend(["--curated-only", "--curation-group", settings.CURATION_GROUP])
        xmlout = self._runNrpsDesigner(args, xmlstr)
        designerDom = parseString(xmlout)
        nrpsList = designerDom.getElementsByTagName('s:component')
        i = 1

        for nrps in nrpsList:
            description = self.description
            if description is not None and len(description) > 0:
                description += "\n"
            elif description is None:
                description = ""
            varname = Substrate.objects.get(pk=monomers[varpos][i]).name
            nrp = NRP.objects.create(owner = self.owner, name="%s %s" % (self.name, varname), description=description + "library variant %d: %s" % (i, varname), indigoidineTagged=self.indigoidineTagged, parent=self)
            for j, monomerId in enumerate(monomers):
                index = 0
                if j == varpos:
                    index = i
                monomer = Substrate.objects.get(pk=monomerId[index])
                SubstrateOrder.objects.create(nrp=nrp, substrate=monomer, order=j)
            nrp._design(self._parseNrps(nrps))
            self.adjustConstruct(nrp)
            i += 1

    def _design(self, domainIdList):
        old_domains = DomainOrder.objects.filter(nrp = self)
        old_domains.delete()
        for count, domainId in enumerate(domainIdList):
            domain = Domain.objects.get(pk = domainId)
            domain_order = DomainOrder.objects.create(nrp= self, domain=domain,
                order=count, left_boundary=None, right_boundary=None)
        self.designed = True
        self.save()
        self.makeConstruct()

    def _parseNrps(self, nrps):
        domainDomList = nrps.getElementsByTagName('s:DnaComponent')
        domainIdList = []
        for comp in domainDomList:
            cid = comp.attributes.getNamedItem("rdf:about").value[1:]
            if cid[:cid.find("_")].isdigit():
                did = comp.getElementsByTagName('s:displayId')[0].firstChild.data
                if did[0] != '_':
                    raise Exception("Could not find domain ID, got %s instead." % did)
                domainIdList.append(int(did[1:]))
        return domainIdList

    def _runNrpsDesigner(self, params, stdin=None):
        lasterror = [None] # nonlocal only in python3
        lastcat = [0]
        logger = logging.getLogger('user_visible')
        xmlout = [""] # nonlocal only in python3
        def processErr(lines):
            def log():
                if lasterror[0]is None:
                    return
                if lastcat[0] == 0:
                    logger.error(lasterror[0])
                elif lastcat[0] == 1:
                    logger.warning(lasterror[0])
                elif lastcat[0] == 2:
                    logger.info(lasterror[0])
            lines = lines.splitlines()
            for line in lines:
                if line.startswith("WARNING: "):
                    log()
                    lastcat[0] = 1
                    lasterror[0] = line[9:]
                elif line.startswith("ERROR: "):
                    log()
                    lastcat[0] = 0
                    lasterror[0] = line[7:]
                elif line.startswith("INFO: "):
                    log()
                    lastcat[0] = 2
                    lasterror[0] = line[6:]
                elif line.startswith(" "):
                    lasterror[0] += line.strip()
                else:
                    log()
                    lastcat[0] = 0
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
        args = ["nrpsdesigner"]
        args.extend(params)
        for k,v in dbSettings.items():
            if v:
                args.extend([k, v])
        child = subprocess.Popen(args, 0, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)

        fds = (child.stdout.fileno(), child.stderr.fileno())
        mask = fcntl.fcntl(fds[0], fcntl.F_GETFL)
        fcntl.fcntl(fds[0], fcntl.F_SETFL, mask | os.O_NONBLOCK)
        mask = fcntl.fcntl(fds[1], fcntl.F_GETFL)
        fcntl.fcntl(fds[1], fcntl.F_SETFL, mask | os.O_NONBLOCK)
        poll = select.poll()
        mask = select.POLLIN | select. POLLPRI | select.POLLERR | select.POLLHUP
        poll.register(fds[0], mask)
        poll.register(fds[1], mask)
        if stdin is not None:
            child.stdin.write(stdin)
            child.stdin.close()
        def processFds(fdsready):
            for fd in fdsready:
                if fd[0] == fds[1] and (fd[1] & select.POLLIN or fd[1] & select.POLLPRI):
                    readFrom(child.stderr, processErr)
                elif fd[0] == fds[0] and (fd[1] & select.POLLIN or fd[1] & select.POLLPRI):
                    readFrom(child.stdout, processOut)
        while child.poll() is None:
            processFds(poll.poll())
        processFds(poll.poll())
        if child.returncode == -11:
            lasterror[0] = "Segmentation fault"
        if child.returncode != 0:
            logger.error("NRPSDesigner returned errorcode %d" % child.returncode)
            raise Exception("NRPSDesigner returned errorcode %d: %s" % (child.returncode, lasterror[0]))
        return xmlout[0]

    def makeSample(self):
        for count, monomerId in enumerate([6,19,14]):
            monomer = Substrate.objects.get(pk=int(monomerId))
            so = SubstrateOrder.objects.create(nrp= self, substrate = monomer, order = count)

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
    left_boundary = models.PositiveIntegerField(null=True, default=None)
    right_boundary = models.PositiveIntegerField(null=True, default=None)

    def __unicode__(self):
        return self.domain.short_name()

