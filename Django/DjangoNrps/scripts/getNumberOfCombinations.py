#!/usr/bin/env python

import os
import sys
import json

if __name__ == "__main__":
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "DjangoNrps.settings")
    if len(sys.argv) < 2:
        print "Need number of monomers"
        exit(0)

    nmonomers = int(sys.argv[1])
    sys.path.append(os.path.dirname(sys.argv[0]) + "/../")
    from designerGui.views import get_available_substrates, toBool
    from databaseInput.models import Substrate
    from django.http import HttpRequest

    if len(sys.argv) > 2:
        curatedonly = toBool(sys.argv[2])
    else:
        curatedonly = True

    combinations = []
    monomers = []
    current = []
    for i in xrange(nmonomers):
        monomers.append([])
        current.append(0)

    monomers[0] = filter(lambda x: x.can_be_added(curatedonly=curatedonly), Substrate.objects.exclude(user__username='sbspks'))

    monomerslist = []
    for i in xrange(1, nmonomers):
        monomerslist.append(monomers[i - 1][current[i - 1]].pk)
        monomers[i] = get_available_substrates(monomerslist, False, curatedonly)[0]

    while current[0] < len(monomers[0]) - 1:
        currcombi = list(monomerslist)
        currcombi.append(monomers[-1][current[-1]])
        combinations.append(currcombi)
        current[-1] += 1

        i = len(current) - 1
        while current[i] >= len(monomers[i]) and i > 0:
            current[i] = 0
            i -= 1
            current[i] += 1
        if i < len(current) - 1:
            monomerslist = [monomers[i][current[i]].pk for i in xrange(nmonomers)]
    print len(combinations)
