# Define here the models for your scraped items
#
# See documentation in:
# http://doc.scrapy.org/topics/items.html

from scrapy.item import Item, Field

class DomainItem(Item):
    sequence = Field()
    dtype = Field()
    substrate = Field()
    start = Field()
    stop = Field()
    linkerbeforestart = Field()
    linkerafterstop = Field()
    module = Field()

class SbspksItem(Item):
    name = Field()
    domains = Field()
    sequences = Field()