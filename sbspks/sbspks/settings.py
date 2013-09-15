# Scrapy settings for sbspks project
#
# For simplicity, this file contains only the most important settings by
# default. All the other settings are documented here:
#
#     http://doc.scrapy.org/topics/settings.html
#

BOT_NAME = 'sbspks'

SPIDER_MODULES = ['sbspks.spiders']
NEWSPIDER_MODULE = 'sbspks.spiders'

DUPEFILTER_CLASS = 'scrapy.dupefilter.BaseDupeFilter'
CONCURRENT_REQUESTS = 200
CONCURRENT_REQUESTS_PER_DOMAIN = 200

ITEM_PIPELINES = ['sbspks.pipelines.SbspksPipeline']

# Crawl responsibly by identifying yourself (and your website) on the user-agent
#USER_AGENT = 'sbspks (+http://www.yourdomain.com)'
