from scrapy.spider import BaseSpider
from scrapy.selector import HtmlXPathSelector
from scrapy.utils.response import get_base_url
from scrapy.utils.url import urljoin_rfc
from scrapy.http import Request

import re

from sbspks.items import SbspksItem, DomainItem
import pdb
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
        for nrps in nrpss:
            link = nrps.select('./@href').extract()[0]
            text = nrps.select('.//text()').extract()[0]
            if text != 'CHLORAMPHENOCOL' and text != 'TYROCIDIN':
                yield Request(urljoin_rfc(base_url, link), callback=self.parseNrps, meta={'name': text.title()})

    def parseNrps(self, response):
        hxs = HtmlXPathSelector(response)
        base_url = get_base_url(response)
        proteins = hxs.select('//a[./text()="FASTA"]')
        module = 1
        item = SbspksItem()
        item['sequences'] = [None] * len(proteins)
        item['sequenceNames'] = [None] * len(proteins)
        item['name'] = response.meta['name']
        item['domains'] = []
        for i in xrange(len(proteins)):
            seq_url = proteins[i].select('./@href').extract()[0]
            yield Request(urljoin_rfc(base_url, seq_url), callback=self.parseSequenceWindow, meta={'sbspksitem': item, 'seqnr': i, 'callback': self.parseProtein})
            if i < len(proteins) - 1:
                domains = proteins[i].select('./following-sibling::a[./text()="FASTA"][1]/preceding-sibling::a[preceding-sibling::a[@href="%s"]]' % seq_url)
            else:
                domains = proteins[i].select('./following-sibling::a')
            ditem = None
            linker = None
            for domain in domains:
                link =  urljoin_rfc(base_url, domain.select('./@href').extract()[0])
                img = domain.select('./img/@src').extract()[0]
                if (img.rfind("linker.gif")) != -1:
                    linker = link
                    if ditem is not None:
                        yield Request(linker, callback=self.parseSequenceWindow, meta={'domain': ditem, 'callback': self.parseLinkerAfter})
                    ditem = None
                else:
                    ditem = DomainItem()
                    ditem['protnr'] = i
                    ditem['dtype'] = re.match("^([A-Za-z]+)", img).group(1).title()
                    if ditem['dtype'] == "C":
                        module += 1
                    elif ditem['dtype'] == "A":
                        substr = domain.select('./following-sibling::text()[1]').extract()[0].strip()
                        if substr[-2:] == "-D":
                            substr = substr[:-2]
                        if substr.lower() in self.names:
                            substr = self.names[substr.lower()]
                        ditem['substrate'] = substr
                        if len(item['domains']) > 0:
                            item['domains'][-1]['substrate'] = substr # for C domain
                    ditem['module'] = module
                    yield Request(link, callback=self.parseSequenceWindow, meta={'domain': ditem, 'callback': self.parseDomain})
                    if linker is not None:
                        yield Request(linker, callback=self.parseSequenceWindow, meta={'domain': ditem, 'callback': self.parseLinkerBefore})
                    linker = None
                    item['domains'].append(ditem)
        yield item

    def parseSequenceWindow(self, response):
        base_url = get_base_url(response)
        hxs = HtmlXPathSelector(response)
        seq = re.search(r"window.open\('([^']+)'", hxs.select('//input[@type="button" and @value="Get Fasta File"]/@onclick').extract()[0]).group(1)
        link = urljoin_rfc(base_url, seq)
        yield Request(link, callback=response.meta['callback'], meta=response.meta)

    def parseProtein(self, response):
        sitem = response.meta['sbspksitem']
        lines = response.body.splitlines()
        sitem['sequences'][response.meta['seqnr']] = self.getSeq(response)
        sitem['sequenceNames'][response.meta['seqnr']] = lines[0]

    def parseDomain(self, response):
        ditem = response.meta['domain']
        ditem['sequence'] = self.getSeq(response)

    def parseLinkerBefore(self, response):
        ditem = response.meta['domain']
        ditem['linkerbeforeseq'] = self.getSeq(response)

    def parseLinkerAfter(self, response):
        ditem = response.meta['domain']
        ditem['linkerafterseq'] = self.getSeq(response)

    def getSeq(self, response):
        lines = response.body.splitlines()
        seq = "".join(lines[1:])
        return re.sub(r'\W+', '', seq) # see e.g. http://www.nii.ac.in/~pksdb/DBASE/webpages/SEQ/bacil002_TE_END_008.seq
