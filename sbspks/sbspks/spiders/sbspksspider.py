from scrapy.spider import BaseSpider
from scrapy.selector import HtmlXPathSelector
from scrapy.utils.response import get_base_url
from scrapy.utils.url import urljoin_rfc

import re

from sbspks.items import SbspksItem, DomainItem

class SBSPKSSpider(BaseSpider):
    name = "sbspks"
    start_urls = [
        "http://www.nii.ac.in/~pksdb/DBASE/listALL4.html",
        ]

    names = {'ala': 'alanine', 'gly': 'Glycine', 'ser': 'serine', 'phe': 'phenylalanine', 'leu': 'leucine', 'thr': 'threonine', 'tyr': 'tyrosine', 'gln': 'glutamine', 'asn': 'asparagine', 'orn': 'ornithine', 'asp': 'aspartic acid', 'glu': 'glutamic acid', 'val': 'valine', 'ile': 'isoleucine', 'trp': 'tryptophane', 'his': 'histidine', 'lys': 'lysine', 'arg': 'arginine', 'met': 'methionine', 'cys': 'cysteine', 'pro': 'proline'}

    def parse(self, response):
        hxs = HtmlXPathSelector(response)
        nrpss = hxs.select('//h1[./font/text()="NRPS"]/following-sibling::p[1]/a')
        base_url = get_base_url(response)
        for nrps in nprss:
            link = nrps.select('./@href').extract()[0]
            text = nrps.select('.//text()').extract()[0]
            if text != 'CHLORAMPHENOCOL' and text != 'TYROCIDIN':
                yield Request(urljoin_rfc(base_url, link), callback=self.parseNrps, meta={'name': text.title()})

    def parseNrps(self, response):
        base_url = get_base_url(response)
        proteins = hxs.select('//a[./text()="FASTA"]')
        module = 1
        item = SbspksItem()
        item['sequences'] = len(proteins)
        item['name'] = response.meta['name']
        item['domains'] = []
        for i in xrange(len(proteins)):
            seq_url = proteins.select('./@href').extract()[0]
            if i < len(proteins) - 1:
                domains = proteins.select('./following-sibling::a[./text()="FASTA"][1]/preceding-sibling::a[preceding-sibling::a[@href="%s"]]' % seq_url)
            else:
                domains = proteins.select('./following-sibling::a')
            ditem = None
            linker = None
            for domain in domains:
                link =  urljoin_rfc(base_url, domain.select('./@href').extract()[0])
                img = domain.select('./img/@src').extract()[0]
                if (img.rfind("linker.gif")) != -1:
                    linker = link
                    if ditem is not None:
                        yield Request(linker, callback=self.parseLinkerAfter, meta={'domain': ditem})
                    ditem = None
                else:
                    ditem = DomainItem()
                    ditem['sequence'] = i
                    ditem['dtype'] = img[:img.find("_")].title()
                    if ditem['dtype'] == "C":
                        module += 1
                    elif ditem['dtype'] == "A":
                        substr = domain.select('./following-sibling::*[1]/text()').extract()[0]
                        if substr[-2:] == "-D":
                            substr = substr[:-2]
                        if substr.lower() in names:
                            substr = names[substr.lower()]
                        ditem['substrate'] = substr
                    ditem['module'] = module
                    yield Request(link, callback=self.parseDomain, meta={'domain': ditem})
                    if linker is not None:
                        yield Request(linker, callback=self.parseLinkerBefore, meta={'domain': ditem})
                    linker = None
                    item['domains'].append(ditem)
        yield item

    def parseDomain(self, response):
        ditem = response.meta['domain']
        limits = self.parseFasta(response)
        ditem['start'] = limits[0]
        ditem['stop'] = limits[1]

    def parseLinkerBefore(self, response):
        ditem = response.meta['domain']
        limits = self.parseFasta(response)
        ditem['linkerbeforestart'] = limits[0]

    def parseLinkerAfter(self, response):
        ditem = response.meta['domain']
        limits = self.parseFasta(response)
        ditem['linkerafterstop'] = limits[1]

    def parseFasta(self, response):
        faheader = response.body.splitlines()[0]
        match = re.search("\\((\d+)\s+(\d+)\\)$", faheader)
        return (int(match.group(1)) * 3, int(match.group(2)) * 3)
