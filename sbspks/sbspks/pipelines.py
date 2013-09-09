# Define your item pipelines here
#
# Don't forget to add your pipeline to the ITEM_PIPELINES setting
# See: http://doc.scrapy.org/topics/item-pipeline.html

from sbspks.items import SbspksItem, DomainItem
from databaseInput.models import *

class SbspksPipeline(object):
    spiders = 0
    data = []

    @classmethod
    def open_spider(cls, spider):
        cls.spiders += 1

    @classmethod
    def close_spider(cls, spider):
        cls.spiders -= 1
        if cls.spiders == 0:
            for item in data:
                prod = Product(name=item['name'])
                prod.save()
                seqs = []
                for i in xrange(item['sequences']):
                    cds = Cds(product=prod)
                    cds.save()
                    seqs.append(cds)
                for i in xrange(len(item['domains'])):
                    d = item['domains'][i]
                    domain = Domain(module=d['module'], cds=seqs[d['sequence']], domainType=d['dtype'], pfamLinkerStart=d['linkerbeforestart'], pfamLinkerStop=d['linkerafterstop'], pfamStart=d['start'], pfamStop=d['stop'])
                    if d['dtype'] == 'A' or d['dtype'] == 'C':
                        try:
                            s = Substrate.objects.get(name="L-%s" % d['substrate'])
                        except Substrate.DoesNotExist:
                            s = Substrate(name=d['substrate'], chirality='L')
                            s.save()
                        domain.substrateSpecificity.add(s)
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
