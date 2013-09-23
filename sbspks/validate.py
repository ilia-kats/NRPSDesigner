#!/usr/bin/env python
# -*- encoding: utf-8 -*-

import sys
import os
import shutil
import random
import inspect

import databaseInput.models
from databaseInput.models import *
from databaseInput.validators import *
from django.contrib.auth.models import User
from django.core.mail import EmailMessage
from django.conf import settings

files = {"createdby" : "createdbysbspks.txt",
         "seqsdefaultorigin" : "seqsdefaultorigin.txt",
         "seqsnoproduct" : "seqsnoproduct.txt",
         "seqsnolinkout" : "seqsnolinkout.txt",
         "seqsinvalid" : "seqsinvalid.txt",
         "deadorigins" : "deadorigins.txt",
         "deadproducts" : "deadproducts.txt",
         "domainsnoseq" : "domainsnoseq.txt",
         "domainsnosubstrate" : "domainsnosubstrate.txt"
        }

def write(fileobj, obj, prepend="", append=""):
    fileobj.write("%sID: %d; %s; %s\n" % (prepend, obj.pk, obj, append))

outdir = "/tmp/%s_%f" % ("nrpsdesignervalidate", random.random())
os.mkdir(outdir)
os.chdir(outdir)

sbspksuser = User.objects.get(username="sbspks")

createdby = open(files["createdby"], "w")
for name, obj in inspect.getmembers(databaseInput.models):
    if hasattr(obj, "user"):
        createdby.write("%s\n" % name)
        for entry in obj.objects.filter(user=sbspksuser):
            write(createdby, entry, "\t")
createdby.close()

origins = Origin.objects.all()
products = Product.objects.all()
defaultorigin = Origin.objects.get(pk=3)
seenorigins = set()
seenproducts = set()
cdss = Cds.objects.all()
seqsdefaultorigin = open(files["seqsdefaultorigin"], "w")
seqsnoproduct = open(files["seqsnoproduct"], "w")
seqsnolinkout = open(files["seqsnolinkout"], "w")
seqsinvalid = open(files["seqsinvalid"], "w")
for cds in cdss:
    if cds.origin == defaultorigin:
        write(seqsdefaultorigin, cds)
    seenorigins.add(cds.origin)
    if cds.product is None or cds.product.name == "":
        write(seqsnoproduct, cds)
    seenproducts.add(cds.product)
    if cds.linkout.count() == 0:
        write(seqsnolinkout, cds)
    try:
        validateCodingSeq(cds.dnaSequence)
    except BaseException as e:
        write(seqsinvalid, cds, append=str(e))
seqsdefaultorigin.close()
seqsnoproduct.close()
seqsinvalid.close()

deadorigins = open(files["deadorigins"], "w")
for ori in origins:
    if ori not in seenorigins:
        write(deadorigins, ori)
deadorigins.close()
deadproducts = open(files["deadproducts"], "w")
for prod in products:
    if prod not in seenproducts:
        write(deadproducts, prod)
deadproducts.close()

dtypeswithsubstr = Type.objects.filter(name__in=['A', 'C'])
domainsnoseq = open(files["domainsnoseq"], "w")
domainsnosubstrate = open(files["domainsnosubstrate"], "w")
for domain in Domain.objects.filter(domainType__in=dtypeswithsubstr):
    if domain.cds is None:
        write(domainsnoseq, domain)
    if domain.substrateSpecificity is None or domain.substrateSpecificity.count() == 0:
        write(domainsnosubstrate, domain)
domainsnoseq.close()
domainsnosubstrate.close()

if len(sys.argv) > 1:
    mail = EmailMessage("NRPSDesigner database validation results", "NRPSDesigner database validation", settings.DEFAULT_FROM_EMAIL, sys.argv[1:])
    for outf in files.values():
        mail.attach_file(outf)
    mail.send()
    shutil.rmtree(outdir)
