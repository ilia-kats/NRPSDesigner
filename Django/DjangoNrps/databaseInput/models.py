from __future__ import unicode_literals

from django.db import models
from django.contrib.contenttypes.models import ContentType
from django.contrib.contenttypes import generic

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

from nrpsSMASH.analyzeNrpCds import nrpsSmash
from databaseInput.MSA.MSA import msa_run
from .validators import validateCodingSeq

class Cds(models.Model):
    origin = models.ForeignKey('Origin')
    product = models.ForeignKey('Product', blank=True, null=True)
    geneName = models.CharField(max_length=100) 
    dnaSequence = models.TextField(validators=[validateCodingSeq])
    description = models.TextField(blank=True, null=True)
    linkout =  generic.GenericRelation('Linkout')
    user = models.ForeignKey('auth.User', blank=True, null=True)
   
    def __unicode__(self):
        return str(self.origin) + self.geneName

    class Meta:
        verbose_name = "Coding sequence"
        verbose_name_plural = "Coding sequences"

    # uses nrpsSMASH stuff in order to generate initial input for formset  
    def predictDomains(self):
        initialDicts = []
        allDomainTypes = [x.smashName for x in Type.objects.all()]
        nrpsSmashResult = nrpsSmash(self.dnaSequence)
        for predictedDomain in nrpsSmashResult.domaindict2['gene']:
            if predictedDomain[0] in allDomainTypes:
                domainDict = {}
                domainDict['pfamStart'] = predictedDomain[1]
                domainDict['pfamStop']  = predictedDomain[2]
                domainType = Type.objects.get(smashName = predictedDomain[0])
                domainDict['domainType'] = domainType
                domainDict['user'] = self.user #possibly improve this afterwards..
                initialDicts.append(domainDict)
        return initialDicts


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
    cds = models.ForeignKey('Cds')
    domainType = models.ForeignKey('Type')
    substrateSpecificity = models.ManyToManyField('Substrate', blank=True, null=True)
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
    linkout =  generic.GenericRelation('Linkout')
    user = models.ForeignKey('auth.User', blank=True, null=True)


    def __unicode__(self):
        return str(self.cds) + str(self.module) + str(self.domainType)

    def get_sequence(self):
        cdsSequence = self.cds.dnaSequence
        domainStart = self.pfamStart
        domainStop  = self.pfamStop
        domainSequence = cdsSequence[domainStart:domainStop]
        return domainSequence

    def get_seq_object(self):
        domainSequence = self.get_sequence()
        domainSeqObject = Seq(domainSequence, IUPAC.unambiguous_dna)
        return domainSeqObject

    def get_seqrecord_object(self):
        domainSeqObject = self.get_seq_object()
        name = self.cds.geneName  + str(self.module) + str(self.domainType)
        domainSeqRecord = SeqRecord(seq=domainSeqObject, name= name, id=name)
        return domainSeqRecord

    def align_same_type(self):
        return self.domainType.align_same_type()


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
        domains = Domain.objects.filter(domainType = self)
        # just b/c db has some empty entries right now
        domains = [dom for dom in domains if len(dom.get_sequence())>0]
        domainSeqRecordList = [dom.get_seqrecord_object() for dom in domains]
        return msa_run(domainSeqRecordList)

class Linkout(models.Model):
    linkoutType = models.ForeignKey('LinkoutType')
    identifier = models.CharField(max_length=50)
    limit = models.Q(app_label = 'databaseInput', model = 'Substrate') | models.Q(app_label = 'databaseInput', model = 'Domain') | models.Q(app_label = 'databaseInput', model = 'Origin') | models.Q(app_label = 'databaseInput', model = 'Cds') | models.Q(app_label = 'databaseInput', model = 'Product')
    content_type = models.ForeignKey(ContentType, limit_choices_to = limit)
    object_id = models.PositiveIntegerField()
    content_object = generic.GenericForeignKey('content_type', 'object_id')
    user = models.ForeignKey('auth.User', blank=True, null=True)

    def __unicode__(self):
        return self.identifier

class LinkoutType(models.Model):
    shortcut = models.CharField(max_length=10)
    url = models.CharField(max_length=100)
    description = models.TextField(blank=True, null=True)

    def __unicode__(self):
        return self.shortcut
