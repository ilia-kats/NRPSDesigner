from __future__ import unicode_literals

from django.db import models
from django.core.exceptions import ObjectDoesNotExist
from django.contrib.contenttypes.models import ContentType
from django.contrib.contenttypes import generic
from django.core.validators import MinLengthValidator
from celery.contrib.methods import task

from itertools import count

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from nrpsSMASH.analyzeNrpCds import nrpsSmash
from databaseInput.MSA.MSA import msa_run
from databaseInput.validators import validateCodingSeq

import logging
import pdb
from time import sleep
class Cds(models.Model):
    origin = models.ForeignKey('Origin')
    product = models.ForeignKey('Product', blank=True, null=True, related_name = "cdsSequence")
    geneName = models.CharField(max_length=100)
    dnaSequence = models.TextField(validators=[validateCodingSeq])
    description = models.TextField(blank=True, null=True)
    linkout =  generic.GenericRelation('Linkout')
    user = models.ForeignKey('auth.User', blank=True, null=True)

    def __unicode__(self):
        return str(self.origin) + " " + self.geneName

    class Meta:
        verbose_name = "Coding sequence"
        verbose_name_plural = "Coding sequences"

    # uses nrpsSMASH stuff in order to generate initial input for formset
    @task()
    def predictDomains(self):
        logging.getLogger('user_visible').info("Automated domain prediction started..")
        initialDicts = []
        allDomainTypes = [x.smashName for x in Type.objects.all()]
        nrpsSmashResult = nrpsSmash(self.dnaSequence)
        try:
            consensusPreds = nrpsSmashResult.consensuspreds
        except AttributeError:
            raise AttributeError("""According to the in-built domain prediction pipeline,
                the DNA sequence you entered does not contain any NRPS domains.""")

        #code below still not nice, should probably adapt nrpsSMASH to give nicer output
        consensusKeys = sorted(consensusPreds, key = lambda x: int(x.split('_A')[-1]))
        consensusValues = [consensusPreds[key] for key in consensusKeys]
        consensusValues = (x for x in consensusValues)
        lasttdomain = None
        for predictedDomain in nrpsSmashResult.domaindict2['gene']:
            if predictedDomain[0] in allDomainTypes:
                domainDict = {}
                domainDict['pfamStart'] = 3*predictedDomain[1]+1
                domainDict['pfamStop']  = 3*predictedDomain[2]+3
                domainType = Type.objects.get(smashName = predictedDomain[0])
                domainDict['domainType'] = domainType
                domainDict['user'] = self.user #possibly improve this afterwards..
                if domainType.name == "A":
                    specificity = consensusValues.next()
                    if specificity in [x.smashName for x in Substrate.objects.all()]:
                        substrate = Substrate.objects.get(smashName = specificity)
                        domainDict['substrateSpecificity'] = [substrate.pk]
                elif domainType.name == "T":
                    lasttdomain = domainDict
                elif (domainType.name == "C_L" or domainType.name == "C_D" or domainType.name == "C") and lasttdomain is not None:
                    lasttdomain['domainType'] = Type.objects.get(name="Tstd")
                    lasttdomain = None
                elif domainType.name == "E" and lasttdomain is not None:
                    lasttdomain['domainType'] = Type.objects.get(name="T_Ep")
                    lasttdomain = None
                initialDicts.append(domainDict)
        module_code = self.type_list_to_modules([x['domainType'].name for x in initialDicts])
        for i,module in enumerate(module_code):
            initialDicts[i]['module'] = module
        return initialDicts

    def get_domain_type_name_list(self):
        return [x.domainType.name for x in self.domains.order_by('pfamStart')]

    # return list of actual domains rather than queryset..
    def get_ordered_domain_list(self):
        return list(self.domains.order_by('pfamStart'))

    @staticmethod
    def type_list_to_modules(type_name_list):
        module_list = []
        c_type_domains = ['C_L','C_D','C_Xg','C_du','C_gl','C_St']
        if type_name_list[0] in c_type_domains:
            counter = count(0)
        else:
            counter = count(1)
        curr_count = counter.next()
        for domainTypeName in (type_name_list):
            if domainTypeName in c_type_domains:
                curr_count = counter.next()
            module_list.append(curr_count)
        return module_list

    # e.g. ATCATTE -> [1,1,2,2,2,2]
    def module_code(self):
        return self.type_list_to_modules(self.get_domain_type_name_list())

    # get DNA sequence based on start and stop domain
    # should be actual domain objects, not IDs
    def get_sequence(self, start_domain, stop_domain, start=None, stop=None):
        if start_domain.cds == self and stop_domain.cds == self:
            domain_list = self.get_ordered_domain_list()
            if start is None or isinstance(start, Domain) or type(start) not in [long, int]:
                startpos = start_domain.get_start(start)
            else:
                startpos = start
            if stop is None or isinstance(stop, Domain) or type(stop) not in [long, int]:
                stoppos = stop_domain.get_stop(stop)
            else:
                stoppos = stop

            return self.dnaSequence[startpos:stoppos]
        else:
            pass  #raise some error..

    def get_biojs_entry(self, protein=False):
        cds_name  = str(self)
        seq   = self.dnaSequence
        if protein:
            seq = translate_dna(seq, cds=True)
        domains = self.domains.all()
        biojs_entry = {'id':cds_name,
            'sequence' : seq,
            'annotations': [x.get_biojs_annotation(protein=protein) for x in domains]
            }
        return biojs_entry


class Origin(models.Model):
    sourceType = models.CharField(max_length=10, choices= (('Species','Species'),('Biobrick','Biobrick'), ('Other', 'Other DNA source')))
    source = models.CharField(max_length=20) #usually taxon ID, can also be biobrick or plasmid identifier
    species = models.CharField(max_length = 100, blank=True, null= True)
    description = models.TextField()
    linkout =  generic.GenericRelation('Linkout')
    parent = models.ForeignKey('self', blank=True, null=True, related_name='child')
    user = models.ForeignKey('auth.User', blank=True, null=True)

    def __unicode__(self):
        if self.sourceType == "Species":
            return self.species
        else:
            return self.source

class Product(models.Model):
    name = models.CharField(max_length=50)
    description = models.TextField()
    linkout = generic.GenericRelation('Linkout')
    user = models.ForeignKey('auth.User', blank=True, null=True)

    def __unicode__(self):
       return self.name

class Domain(models.Model):
    module = models.IntegerField()
    cds = models.ForeignKey('Cds', related_name = "domains")
    domainType = models.ForeignKey('Type')
    substrateSpecificity = models.ManyToManyField('Substrate', blank=True, null=True, related_name = "adenylationDomain")
    chirality = models.CharField(max_length=1, choices= (('L','L'),('D','D'),('N','None')), default='N')
    description = models.TextField(blank=True, null=True)
    pfamLinkerStart = models.IntegerField(blank=True, null=True)
    pfamLinkerStop = models.IntegerField(blank=True, null=True)
    definedLinkerStart = models.IntegerField(blank=True, null=True)
    definedLinkerStop = models.IntegerField(blank=True, null=True)
    pfamStart = models.IntegerField()
    pfamStop = models.IntegerField()
    definedStart = models.IntegerField(blank=True, null=True)
    definedStop = models.IntegerField(blank=True, null=True)
    next_domain = models.ManyToManyField('self', symmetrical=False, through='DomainTuple', blank=True, related_name='prev_domain')
    linkout =  generic.GenericRelation('Linkout')
    user = models.ForeignKey('auth.User', blank=True, null=True)


    def __unicode__(self):
        return str(self.cds) + str(self.module) + str(self.domainType)

    def short_name(self):
        return str(self.cds.geneName) + str(self.module) + str(self.domainType)


    def get_start(self, prevd=None, with_linker=True):
        """returns in python format: 0-based"""
        start = None
        if prevd is not None:
            try:
                # TODO: need to find a better way to handle multiple experimental data
                start = self.prev_domain.get(pk=prevd.pk).next_tuple.filter(next_domain__pk=self.pk)[0].next_position
            except ObjectDoesNotExist:
                pass
        if with_linker:
            if self.definedLinkerStart is not None:
                start = self.definedLinkerStart
            elif self.pfamLinkerStart is not None:
                start = self.pfamLinkerStart
        if self.definedStart is not None:
            start = self.definedStart
        else:
            start = self.pfamStart
        return start - 1

    def get_stop(self, nextd=None, with_linker=False):
        """returns in python format: 0-based, 1 position after end"""
        stop = None
        if nextd is not None:
            try:
                # TODO: need to find a better way to handle multiple experimental data
                stop = self.next_domain.get(pk=nextd.pk).prev_tuple.filter(prev_domain__pk=self.pk)[0].prev_position
            except ObjectDoesNotExist:
                pass
        if with_linker:
            if self.definedLinkerStop is not None:
                stop = self.definedLinkerStop
            elif self.pfamLinkerStop is not None:
                stop = self.pfamLinkerStop
        if self.definedStop is not None:
            stop = self.definedStop
        else:
            stop = self.pfamStop
        return stop

    def get_sequence(self, includeForwardLinker=False, includeAftLinker=False):
        cdsSequence = self.cds.dnaSequence
        domainStart = self.get_start(with_linker=includeForwardLinker)
        domainStop  = self.get_stop(with_linker=includeAftLinker)
        domainSequence = cdsSequence[domainStart:domainStop]
        return domainSequence

    def get_seq_object(self, protein= False):
        domainSequence = self.get_sequence()
        domainSeqObject = Seq(domainSequence, IUPAC.unambiguous_dna)
        if protein:
            domainSeqObject = domainSeqObject.translate()
        return domainSeqObject

    def get_seqrecord_object(self, protein=False):
        domainSeqObject = self.get_seq_object(protein=protein)
        name = self.cds.geneName  + str(self.module) + str(self.domainType)
        domainSeqRecord = SeqRecord(seq=domainSeqObject, name= name, id=name)
        return domainSeqRecord

    @task()
    def align_same_type(self):
        return (self.pk, self.domainType.align_same_type())

    # introduce this function, so we can deduce which A domains are
    # specific for only 1 substrate
    def number_of_specificities(self):
        return len(self.substrateSpecificity.all())

    def get_biojs_annotation(self, protein=False):
        short_name  = self.short_name()
        name = str(self)
        start = self.get_start() + 1
        stop  = self.get_stop()
        if protein:
            start, stop = dna_to_prot_coords(start, stop)
        annotation = {'name': short_name,
            'html': name,
            'regions': [{'start':start, 'end':stop,'color':"blue"}]
            }
        return annotation

    def get_biojs_highlight(self, protein=False):
        start = self.get_start() + 1
        stop  = self.get_stop()

        if protein:
            start, stop = dna_to_prot_coords(start,stop)
        highlight = { 'start':start ,
            'end':stop,
            'color':"white",
            'background':"green"}
        return highlight
    # return python dict to be converted to json for use by bioJs sequence
    def get_biojs_entry(self, protein=False):
        biojs_entry = self.cds.get_biojs_entry(protein=protein)
        biojs_entry['highlights'] = [self.get_biojs_highlight(protein=protein)]
        return biojs_entry


class Substrate(models.Model):
    name = models.CharField(max_length=30)
    chirality = models.CharField(max_length=1, choices= (('L','L'),('D','D')))
    structure = models.TextField()
    linkout =  generic.GenericRelation('Linkout')
    enantiomer = models.ForeignKey('self', blank=True, null=True)
    modification = models.ManyToManyField('Modification', blank=True, null=True)
    parent = models.ForeignKey('self', blank=True, null=True, related_name='child')
    user = models.ForeignKey('auth.User', blank=True, null=True)
    smashName = models.CharField(max_length=50, blank=True, null=True)

    def __unicode__(self):
        return self.name

    def reverse_chirality(self):
        if self.chirality == "L":
            return "D"
        else:
            return "L"

    def can_be_added_by_adenylation_domain(self, curatedonly=False):
        domains = self.adenylationDomain.annotate(models.Count('substrateSpecificity')).exclude(substrateSpecificity__count__gt=1, user__username='sbspks').filter(domainType__name='A')
        if curatedonly:
            domains = domains.filter(user__groups__name='curator')
        return domains.count() > 0

    def can_be_added_by_condensation_adenylation_domain(self, chirality, curatedonly=False):
        domains =  self.adenylationDomain.annotate(models.Count('substrateSpecificity')).exclude(substrateSpecificity__count__gt=1, user__username='sbspks')
        if curatedonly:
            domains = domains.filter(user__groups__name='curator')
        return domains.filter(domainType__name='A').count() > 0 and domains.filter(domainType__name='C_' + chirality).count() > 0

    def can_be_added_by_modification_domain(self):
        bla = [self.parent is not None]
        pass

    def can_be_added(self, previousChirality=None, curatedonly=False):
        if previousChirality is None:
            if self.can_be_added_by_adenylation_domain(curatedonly): #or self.can_be_added_by_modification_domain():
                return True
            elif self.enantiomer is not None:
                return self.enantiomer.can_be_added_by_adenylation_domain(curatedonly) #or self.enantiomer.can_be_added_by_modification_domain()
            else:
                return False
        else:
            if self.can_be_added_by_condensation_adenylation_domain(previousChirality, curatedonly): #or self.can_be_added_by_modification_domain():
                return True
            elif self.enantiomer is not None:
                return self.enantiomer.can_be_added_by_condensation_adenylation_domain(previousChirality, curatedonly) #or self.enantiomer.can_be_added_by_modification_domain()
            else:
                return False

class Modification(models.Model):
    name = models.CharField(max_length=100)
    domainType =  models.ForeignKey('Type', blank=True, null=True, related_name='modificationAdded')

    def __unicode__(self):
        return self.name

class Type(models.Model):
    name = models.CharField(max_length=4)
    isModification = models.BooleanField()
    pfamName = models.CharField(max_length=100, blank=True, null=True)
    pfamId   = models.CharField(max_length=20, blank=True, null=True)
    description = models.TextField()
    pfamGraphic = models.TextField()
    smashName = models.CharField(max_length=50, blank=True, null=True)

    def __unicode__(self):
        return self.name

    def align_same_type(self):
        ttypes = ["T", "Tstd", "T_Ep"]
        if self.name in ttypes:
            types = []
            for ttype in ttypes:
                types.append(Type.objects.get(name=ttype))
            domains = Domain.objects.filter(domainType__in=types)
        else:
            domains = Domain.objects.filter(domainType = self)
        # just b/c db has some empty entries right now
        domains = [dom for dom in domains if len(dom.get_sequence())>0]
        domainSeqRecordList = [dom.get_seqrecord_object(protein=True) for dom in domains]
        return msa_run(domainSeqRecordList)

# models experimentally validated domain combinations
class DomainTuple(models.Model):
    prev_domain    = models.ForeignKey('Domain', related_name = "next_tuple",verbose_name="left domain")
    next_domain    = models.ForeignKey('Domain', related_name = "prev_tuple",verbose_name="right domain")
    prev_position  = models.IntegerField()
    next_position  = models.IntegerField()
    experiment     = models.ForeignKey('Experiment' , related_name="domain_tuples")

    def __unicode__(self):
        return str(self.prev_domain)+ "with" + str(self.next_domain)

class Experiment(models.Model):
    description = models.TextField()
    linkout =  generic.GenericRelation('Linkout')
    user= models.ForeignKey('auth.User', blank=True, null=True)

    def __unicode__(self):
        return self.description[0:20]



class LinkoutType(models.Model):
    shortcut = models.CharField(max_length=10)
    url = models.CharField(max_length=100)
    description = models.TextField(blank=True, null=True)

    def __unicode__(self):
        return self.shortcut

def get_pubmed_type():
    return LinkoutType.objects.get(shortcut="PubMed");


class Linkout(models.Model):
    linkoutType = models.ForeignKey('LinkoutType')
    identifier = models.CharField(max_length=50, validators=[MinLengthValidator(3)])
    limit = models.Q(app_label = 'databaseInput', model = 'Substrate') | models.Q(app_label = 'databaseInput', model = 'Domain') | models.Q(app_label = 'databaseInput', model = 'Origin') | models.Q(app_label = 'databaseInput', model = 'Cds') | models.Q(app_label = 'databaseInput', model = 'Product') | models.Q(app_label = 'databaseInput', model = 'Experiment')
    content_type = models.ForeignKey(ContentType, limit_choices_to = limit)
    object_id = models.PositiveIntegerField()
    content_object = generic.GenericForeignKey('content_type', 'object_id')
    user = models.ForeignKey('auth.User', blank=True, null=True)

    def __unicode__(self):
        return self.identifier



# some dna-protein helper functions

def translate_dna(dna_seq, cds=True):
    prot_seq = str(Seq(dna_seq, IUPAC.unambiguous_dna).translate(cds=cds, table="Bacterial"))
    return prot_seq

def dna_to_prot_coords(start,stop):
    if start % 3 == 1 and stop % 3 ==0:
        prot_start = (start+2)/3
        prot_stop = stop/3
        return (prot_start,prot_stop)
    else:
        raise ValueError("Can't convert DNA coords to protein coords")

def prot_to_dna_coords(start,stop):
    return (3*start-2, 3*stop)
