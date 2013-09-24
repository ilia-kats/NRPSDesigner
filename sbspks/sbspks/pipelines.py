# Define your item pipelines here
#
# Don't forget to add your pipeline to the ITEM_PIPELINES setting
# See: http://doc.scrapy.org/topics/item-pipeline.html

from sbspks.items import SbspksItem, DomainItem
from databaseInput.models import *
from django.contrib.auth.models import User

import pdb
class SbspksPipeline(object):
    spiders = 0
    data = []

    names = {'ala': 'alanine', 'gly': 'Glycine', 'ser': 'serine', 'phe': 'phenylalanine', 'leu': 'leucine', 'thr': 'threonine', 'tyr': 'tyrosine', 'gln': 'glutamine', 'asn': 'asparagine', 'orn': 'ornithine', 'asp': 'aspartic acid', 'glu': 'glutamic acid', 'val': 'valine', 'ile': 'isoleucine', 'trp': 'tryptophan', 'his': 'histidine', 'lys': 'lysine', 'arg': 'arginine', 'met': 'methionine', 'cys': 'cysteine', 'pro': 'proline'}

    @classmethod
    def open_spider(cls, spider):
        cls.spiders += 1

    @classmethod
    def close_spider(cls, spider):
        cls.spiders -= 1
        if cls.spiders == 0:
            dummyori = Origin.objects.get(pk=3)
            user = User.objects.get(username='sbspks')
            types = {}
            for item in cls.data:
                prod = Product(name=item['name'], user=user)
                prod.save()
                seqs = []
                for i in xrange(len(item['sequences'])):
                    cds = Cds(product=prod, origin=dummyori, geneName=item['sequenceNames'][i][:100], user=user)
                    cds.save()
                    seqs.append(cds)
                for i in xrange(len(item['domains'])):
                    d = item['domains'][i]
                    seq = item['sequences'][d['protnr']]
                    try:
                        linkerbeforestart = seq.index(d['linkerbeforeseq']) + 1
                        linkerbeforestop = linkerbeforestart + len(d['linkerbeforeseq'])
                        domainstart = seq.index(d['sequence']) + 1
                        domainstop = domainstart + len(d['sequence'])
                        linkerafterstart = seq.index(d['linkerafterseq']) + 1
                        linkerafterstop = linkerafterstart + len(d['linkerafterseq'])
                    except ValueError as e:
                        print item
                        print d
                        raise e
                    try:
                        try:
                            dtype = Type.objects.get(name=d['dtype'])
                        except Type.DoesNotExist:
                            dtype = Type(name=d['dtype'])
                            dtype.save()
                        domain = Domain(module=d['module'], cds=seqs[d['protnr']], domainType=dtype, pfamLinkerStart=linkerbeforestart, pfamLinkerStop=linkerafterstop, pfamStart=domainstart, pfamStop=domainstop, definedLinkerStart=linkerbeforestart, definedLinkerStop=linkerafterstop, definedStart=domainstart, definedStop=domainstop, user=user)
                    except BaseException as e:
                        print d
                        print item
                        raise e
                    domain.save() # need pk for many-to-many
                    if d['dtype'] == 'A' or d['dtype'] == 'C':
                        try:
                            if d['substrate'] in cls.names:
                                if d['substrate'] != 'gly':
                                    substr = "L-%s" % cls.names[d['substrate']]
                                else:
                                    substr = cls.names[d['substrate']]
                                s = Substrate.objects.get(name=substr)
                            else:
                                substr = d['substrate']
                                s = Substrate(name='L-%s' % d['substrate'], chirality='L', user=user)
                                s.save()
                                sd = Substrate(name='D-%s' %d['substrate'], chirality='D', user=user)
                                sd.save()
                                s.enantiomer = sd
                                sd.enantiomer = s
                                s.save()
                                sd.save()
                                cls.names[d['substrate']] = d['substrate']
                                domain.substrateSpecificity.add(s)
                        except KeyError as e:
                            print "KeyError:"
                            print d
                            print item

                    if d['dtype'] == 'C':
                        if item['domains'][i - 1]['dtype'] == 'E':
                            domain.chirality = 'D'
                        else:
                            domain.chirality = 'L'
                    domain.save()

    @classmethod
    def process_item(cls, item, spider):
        if isinstance(item, SbspksItem):
            cls.data.append(item)
        return item
