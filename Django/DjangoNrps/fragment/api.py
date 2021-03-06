"""Server-side Fragment Library API"""
#get/set = view
#read/write = database operation

from fragment.models import *
from django.http import HttpResponse, HttpResponseRedirect, HttpResponseNotFound, Http404
from django.shortcuts import render_to_response
from django.contrib.auth.decorators import login_required
from django.contrib.auth.models import Group
from django.views.decorators.http import condition
from django.core.exceptions import ObjectDoesNotExist
from django.db.models import Q
from django.conf import settings

from Bio.Alphabet import IUPAC
from Bio.Seq import Seq

import json
from gibson.jsonresponses import JsonResponse, RawJsonResponse, ERROR

## Helpful functions

the_alphabet = {
        'A': 'T',
        'B': 'V',
        'C': 'G',
        'D': 'H',
        'G': 'C',
        'H': 'D',
        'K': 'M',
        'M': 'K',
        'N': 'N',
        'R': 'Y',
        'S': 'S',
        'T': 'A',
        'V': 'B',
        'W': 'W',
        'Y': 'R',
        'a': 't',
        'b': 'v',
        'c': 'g',
        'd': 'h',
        'g': 'c',
        'h': 'd',
        'k': 'm',
        'm': 'k',
        'n': 'n',
        'r': 'y',
        's': 's',
        't': 'a',
        'v': 'b',
        'w': 'w',
        'y': 'r',
        }

def get_alphabet():
    return the_alphabet

def get_gene(usr, fid):
    if usr.is_authenticated():
        user = usr
    else:
        user = None
    try:
        fid = int(fid)
    except ValueError:
        raise Http404
    try:
        g = Gene.objects.get(Q(owner=user) | Q(viewable='G'), id = fid)
    except ObjectDoesNotExist:
        return False
    return g

def get_gene_info(g):
    info = {
        'id': g.id,
        'name': g.name,
        'desc': g.description,
        'length': g.length(),
        'viewable': g.viewable
    }
    try:
        domaintypes = []
        domains = g.domaingene.domains.all()
        for domain in domains:
            substrates = []
            for substrate in domain.substrateSpecificity.all():
                substrates.append({'name': substrate.name, 'chirality': substrate.chirality})
            jsonlinks = []
            for l in domain.linkout.all():
                jsonlinks.append({'name': l.linkoutType.shortcut, 'url': l.linkoutType.url % l.identifier})
            domaintypes.append({'type': domain.domainType.name, 'substrates': substrates, 'chirality': domain.chirality, 'description': domain.description, 'curated': len(domain.user.groups.filter(name=settings.CURATION_GROUP)) > 0, 'module': domain.module, 'linkouts': jsonlinks})
        info['domaingene'] = True
        info['domaintypes'] = domaintypes
        info['gene'] = domains[0].cds.geneName
        info['species'] = domains[0].cds.origin.species
        info['source'] = domains[0].cds.origin.source
        info['sourceType'] = domains[0].cds.origin.sourceType
    except DomainGene.DoesNotExist:
        pass
    return info

def read_meta(g):
    """Return JSON-ready metadata about a fragment"""
    refs = []
    for r in g.references.all():
        refs.append( {
            'title': r.title,
            'authors': r.authors,
            'journal': r.journal,
            'medline_id': r.medline_id,
            'pubmed_id': r.pubmed_id,
            })
    annots = {}
    for a in g.annotations.all():
        annots[a.key] = a.value
    info = get_gene_info(g)
    info['refs'] = refs
    info['annots'] = annots
    info['origin'] = g.get_origin_display()
    info['length'] = len(g.sequence)
    return info

def write_meta(g, meta):
    """save meta to g"""
    g.name = meta.get('name', g.name)
    g.description = meta.get('desc', g.description)
    #only change references if they've been given
    if meta.has_key('refs'):
        #clear out previous references
        Reference.remove(g)
        Reference.add(g, meta['refs'])
    if meta.has_key('annots'):
        #clear out previous annotations
        Annotation.remove(g)
        for (key, value) in meta['annots']:
            Annotation.add(g, key, value)
    try:
        g.save()
    except Exception as e:
        print "raised exception of type %s: %s" % (type(e), e)

def read_features(g):
    """Read all the features of a fragment"""
    data = []
    for f in g.features.all():
        quals = []
        for q in f.qualifiers.all():
            quals.append({
                'name': q.name,
                'data': q.data,
                })
        s = None
        if f.direction == 'f':
            s = 1
        elif f.direction == 'r':
            s = -1
        data.append({
            'start': f.start,
            'end': f.end,
            'strand': s,
            'type': f.type,
            'id': f.id,
            'qualifiers': quals,
            })
    return data

def write_features(g, feats):
    """Update the features in the database"""
    #remove old features
    Feature.remove(g)
    for f in feats:
        #make sure data makes sense
        start = int(f.get('start',0))
        end = int(f.get('end',0))
        strand = int(f.get('strand', 1))
        if strand not in [-1, 1]:
            strand = 1
        if strand == -1 and start < end:
            t = end
            end = start
            start = t
        elif strand == 1 and start > end:
            t = end
            end = start
            start = t

        quals = {}
        for q in f.get('qualifiers'):
            quals[q.get('name')] = q.get('value')

            ft = SeqFeature(
                    location=FeatureLoaction(start, end),
                    type=f.get('type', ''),
                    strand=strand,
                    id=f.get('id'),
                    qualifiers=quals,
                    )

            Feature.add(ft, g)

def chunk_sequence(g, chunk_size, pad_char):
    """return the sequence in chunks of size chunk_size"""
    i = 0
    length = len(g.sequence)
    while (i+chunk_size) < length:
        yield g.sequence[i:i+chunk_size]
        yield pad_char * 1024
        i = i + chunk_size
    #return the last piece
    yield g.sequence[i:]

def set_part_types(request):
    if request.method == 'POST':
        types = json.loads(request.POST['part_types'])
        for t in types:
            if t['part_type'] != '':
                g = get_gene(request.user, t['id'])
                if g:
                    try:
                        a = g.annotations.get(key='part_type')
                        a.value = t['part_type']
                    except ObjectDoesNotExist:
                        a = Annotation(gene=g, key='part_type', value=t['part_type'])
                    a.save()
        return JsonResponse('OK')

# Actual API stuff
def get_fragment(request, fid):
    """Return a small amount of information about the fragment"""
    try:
        fid = int(fid)
    except ValueError:
        return JsonResponse("Invalid fragment id: %s" % fid, ERROR)
    g = get_gene(request.user, fid)
    if g:
        return JsonResponse(get_gene_info(g))
    else:
        return HttpResponseNotFound()

def list_fragments(request):
    """Return metadata of all fragments owned by a user"""
    if request.user.is_authenticated():
        args = Q(owner=request.user, viewable__in=['L', 'G']) | Q(viewable__in=['G'])
    else:
        args = Q(viewable__in = ['G'])
    try:
        if 'type[]' in request.POST:
            args = args & Q(annotations__key__exact='part_type', annotations__value__in=request.POST.getlist('type[]'))
            frags = Gene.objects.select_related('annotations').filter(args)
        else:
            frags = Gene.objects.filter(args)
    except ObjectDoesNotExist:
        return JsonResponse('Could not read fragments', ERROR)
    ret = []
    for f in frags:
        ret.append(read_meta(f))
    return JsonResponse(ret)

def get_meta(request, fid):
    """get a particular fragment's metadata"""
    g = get_gene(request.user, fid)
    if g:
        return JsonResponse( read_meta(g))
    else:
        return HttpResponseNotFound()

def set_meta(request, fid):
    """Update a fragment's metadata"""
    if request.method == 'POST':
        g = get_gene(request.user, fid)
        if g:
            meta = {'name': request.POST.get('name'),
                            'desc': request.POST.get('desc'),
                            }
            meta['annots'] = []
            for i in range(0, len(request.POST) -2):
                m = request.POST.getlist('annots[%s][]'%i)
                if m:
                    meta['annots'].append(m)
            print 'meta = %s' % meta
            if meta:
                write_meta(g, meta)
                return get_meta(request, fid)
            return JsonResponse("No metadata supplied.", ERROR)
        else:
            return JsonResponse("Fragment with ID='%s' does not exist." % fid, ERROR)
    raise Http404

def get_features(request, fid):
    """Get a fragment's features"""
    g = get_gene(request.user, fid)
    if g:
        return JsonResponse({
            'feats': read_features(g),
            'alpha': get_alphabet(),
            'length': len(g.sequence),
        })
    else:
        return HttpResponseNotFound()

def set_features(request, fid):
    """Save a fragment's features"""
    if request.method == 'POST':
        g = get_gene(request.user, fid)
        if g:
            feats = request.POST.get('features')
            if feats:
                write_features(g, feats)
                return JsonResponse('Done')
            return JsonResponse('No features provided', ERROR)
        else:
            JsonResponse("Fragment with ID='%s' does not exist" % fid, ERROR)
    raise Http404

@condition(etag_func=None)
def get_sequence(request, fid):
    """Stream the sequence in blocks of 1kB by default"""
    try:
        chunk_size = int(request.GET.get('chunk_size', 1024))
    except ValueError:
        return JsonResponse("ERROR: Invalid chunk_size '%s'." %
                request.GET.get('chunk_size', 1024), ERROR)

    pad_char = request.GET.get('pad_char', ' ')
    g = get_gene(request.user, fid)
    if g:
        resp = HttpResponse(chunk_sequence(g, chunk_size, pad_char), mimetype="text/plain")
        return resp
    return JsonResponse('Could Not Stream Sequence', ERROR)

def crop(request, fid):
    """crop the fragment as specified"""
    try:
        fid = int(fid)
    except ValueError:
        return JsonResponse("Invalid fragment id: %s" % fid, ERROR)
    g = get_gene(request.user, fid)
    if not g:
        return JsonResponse('Could not get fragment %s'%fid, ERROR)
    try:
        start = int(request.POST.get('start'))
        end = int(request.POST.get('end'))
    except ValueError:
        return JsonResponse('Start and end must be valid integers', ERROR)

    try:
        f_internal = bool(int(request.POST.get('f_internal', 1)))
        f_all = bool(int(request.POST.get('f_all', 0)))
    except ValueError:
        return JsonResponse('f_internal and f_all must be 1 or 0', ERROR)
    if f_all and not f_internal:
        return JsonResponse(
            'Cannot keep all features without keeping internal ones', ERROR)
    result = request.POST.get('result', 'new')
    if result not in ['new', 'overwrite']:
        return JsonResponse(
            'Result type must be \'new\' or \'overwrite\', not %s' % result,
                ERROR)
    new_name = request.POST.get('new_name')
    new_desc = request.POST.get('new_desc')
    if (result == 'new') and not new_name:
        return JsonResponse(
            'Must specify non-empty new_name when result type is new', ERROR)

    #build a list of features we want to keep
    feats = []
    if f_internal:
        for f in g.features.all():
            #internal features
            if (f.start >= start) and (f.end <= end):
                feats.append(f)
            elif f_all:
                if ( ((f.start > start) and (f.start < end)) or #starts in selection
                        ((f.end < end) and (f.end > start)) or #ends in selection
                        ((f.start < start) and (f.end > end)) ): #selection is subset
                    feats.append(f)

    #get or make the target gene
    if result == 'new':
        target = Gene(owner=request.user,
            name=new_name,
            description=new_desc,
            sequence=g.sequence[start:end])
        target.save()
        #copy references
        for r in g.references.all():
            nr = Reference(gene=target,
                title=r.title,
                authors=r.authors,
                journal=r.journal,
                medline_id=r.medline_id,
                pubmed_id=r.pubmed_id)
            nr.save()
        #copy annotations
        for a in g.annotations.all():
            na = Annotation(gene=target,
                key=a.key,
                value=a.value)
            na.save()
    elif result == 'overwrite':
        target = g
        g.sequence = g.sequence[start:end]
        #clear out old features
        Feature.remove(g)

    #copy the features for adding onto the new fragment
    for f in feats:
        nf = Feature(gene=target,
            type=f.type,
            direction=f.direction,
            start=max(f.start - start, 0),
            end=min(f.end - start, end - start))
        nf.save()

    target.save()

    return JsonResponse(target.id)
