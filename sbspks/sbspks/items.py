# Define here the models for your scraped items
#
# See documentation in:
# http://doc.scrapy.org/topics/items.html

from scrapy.item import Item, Field

class DomainItem(Item):
    sequence = Field()
    protnr = Field()
    dtype = Field()
    substrate = Field()
    linkerbeforeseq = Field()
    linkerafterseq = Field()
    module = Field()

class SbspksItem(Item):
    name = Field()
    domains = Field()
    sequences = Field()
    sequenceNames = Field()
