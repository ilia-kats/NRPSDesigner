from os import path
import os
from argparse import Namespace
from multiprocessing import cpu_count
import subprocess
import logging
import shutil
import tempfile
import inspect
try:
    from cStringIO import StringIO
except ImportError:
    from StringIO import StringIO
# Don't display the SearchIO experimental warning, we know this.
import warnings
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    from Bio import SearchIO

# config example
# Namespace(all_orfs=False, check_prereqs_only=False, 
# clusterblast=False, cpus=8, debug=False,
# enabled_cluster_types=[ 'thiopeptide', 'phosphonate', 'other'], 
# end=-1, eukaryotic=False, full_blast=False, 
# full_hmmer=False, fullhmmer=Namespace(evalue='0.01', score='0'), 
# glimmer=Namespace(basedir='/usr/lib/tigr-glimmer', 
# ignore_score_length='3000', max_overlap='50', min_gene_length='90',
# threshold='30', translation_table='11'), 
# inclusive=False, input_type='nucl', 
# list_plugins=False, 
# sequences=['tycA.gb'], smcogs=False, 
# start=-1, subclusterblast=False, verbose=False)


# hmmscan: Search a protein sequence against a protein profile HMM database.
def run_hmmscan(target_hmmfile, query_sequence, opts=None):
    "Run hmmscan"
    # think about this afterwards..
    #config = get_config()
    config = Namespace()
    config.cpus = cpu_count()
    command = ["hmmscan", "--cpu", str(config.cpus)]
    if opts is not None:
        command.extend(opts)
    command.extend([target_hmmfile, '-'])
    try:
        out, err, retcode = execute(command, input=query_sequence)
    except OSError:
        return []
    if retcode != 0:
        logging.debug('hmmscan returned %d: %r while scanning %r' % (retcode,
                        err, query_sequence))
        return []
    res_stream = StringIO(out)
    results = list(SearchIO.parse(res_stream, 'hmmer3-text'))
    return results

def hmmlengths(hmmfile):
    hmmlengthsdict = {}
    openedhmmfile = open(hmmfile,"r")
    filetext = openedhmmfile.read()
    filetext = filetext.replace("\r","\n")
    hmms = filetext.split("//")[:-1]
    for i in hmms:
        namepart = i.split("NAME  ")[1]
        name = namepart.split("\n")[0]
        lengthpart = i.split("LENG  ")[1]
        length = lengthpart.split("\n")[0]
        hmmlengthsdict[name] = int(length)
    return hmmlengthsdict

def execute(commands, input=None):
    "Execute commands in a system-independent manner"

    if input is not None:
        stdin_redir = subprocess.PIPE
    else:
        stdin_redir = None

    proc = subprocess.Popen(commands, stdin=stdin_redir,
                            stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE)
    try:
        out, err = proc.communicate(input=input)
        retcode = proc.returncode
        return out, err, retcode
    except OSError, e:
        logging.debug("%r %r returned %r" % (commands, input[:40], e))
        raise

### HELPER CLASSES ###

class Storage(dict):
    """Simple storage class"""
    def __init__(self, indict=None):
        if indict is None:
            indict = {}
        dict.__init__(self, indict)
        self.__initialized = True

    def __getattr__(self, attr):
        try:
            return self.__getitem__(attr)
        except KeyError:
            raise AttributeError(attr)

    def __setattr__(self, attr, value):
        if not self.__dict__.has_key('_Storage__initialized'):
            return dict.__setattr__(self, attr, value)
        elif attr in self:
            dict.__setattr__(self, attr, value)
        else:
            self.__setitem__(attr, value)

# from https://github.com/kblin/bioinf-helperlibs/blob/master/helperlibs/wrappers/io.py

class TemporaryDirectory(object):
    def __init__(self, suffix='', prefix='tmp', dir=None, change=False):
        self.change = change
        self.tempdir = tempfile.mkdtemp(suffix, prefix, dir)

    def __enter__(self):
        if self.change:
            self.old_wd = os.getcwd()
            os.chdir(self.tempdir)
        return self.tempdir

    def __exit__(self, type, value, traceback):
        if self.change:
            os.chdir(self.old_wd)
        shutil.rmtree(self.tempdir)


def get_full_path(current_file, file_to_add):
    "Get the full path of file_to_add in the same directory as current_file"
    return path.join(path.dirname(path.abspath(current_file)), file_to_add)

def writefasta(names, seqs, filename):
    "Write sequence to a file"
    e = 0
    f = len(names) - 1
    out_file = open(filename,"w")
    while e <= f:
        out_file.write(">")
        out_file.write(names[e])
        out_file.write("\n")
        out_file.write(seqs[e])
        out_file.write("\n")
        e += 1
    out_file.close()

# LALALLALALALLA
# some functions to work with SeqRecords
# gotta check later which should be removed from our codebase
# LALLALALALALLALA

def get_all_features_of_type(seq_record, types):
    "Return all features of the specified types for a seq_record"
    if isinstance(types, str):
        # force into a tuple
        types = (types, )
    features = []
    for f in seq_record.features:
        if f.type in types:
            features.append(f)
    return features

def get_feature_dict(seq_record):
    """Get a dictionary mapping features to their IDs"""
    features = get_cds_features(seq_record)
    feature_by_id = {}
    for feature in features:
        gene_id = get_gene_id(feature)
        feature_by_id[gene_id] = feature
    return feature_by_id

def get_cds_features(seq_record):
    "Return all CDS features for a seq_record"
    return get_all_features_of_type(seq_record, "CDS")

def get_gene_id(feature):
    "Get the gene ID from locus_tag, gene name or protein id, in that order"
    if 'locus_tag' in feature.qualifiers:
        return feature.qualifiers['locus_tag'][0]
    if 'gene' in feature.qualifiers:
        return feature.qualifiers['gene'][0]
    if 'protein_id' in feature.qualifiers:
        return feature.qualifiers['protein_id'][0]
    return "no_tag_found"

def get_pksnrps_cds_features(seq_record):
    features = get_cds_features(seq_record)
    for feature in features:
        if not 'sec_met' in feature.qualifiers or len([feat for feat in feature.qualifiers['sec_met'] if "Type: " in feat]) == 0:
            continue
    secmetgenes = [feat for feat in features if 'sec_met' in feat.qualifiers]
    pksnrpscoregenes = [feat for feat in secmetgenes if len([qual for qual in feat.qualifiers['sec_met'] if "NRPS/PKS Domain:" in qual]) > 0]
    return pksnrpscoregenes

def get_aa_sequence(seq_record, feature, to_stop=False):
    """Extract sequence from specific CDS feature in sequence record"""
    return str(seq_record.seq.translate(to_stop=True))
    #fasta_seq = feature.qualifiers['translation'][0]
    #if "*" in fasta_seq:
    #     if to_stop:
    #         fasta_seq = fasta_seq.split('*')[0]
    #     else:
    #         fasta_seq = fasta_seq.replace("*","X")
    # if "-" in fasta_seq:
    #     fasta_seq = fasta_seq.replace("-","")
    # return fasta_seq

def get_specific_multifasta(seq_record, features):
    """Extract multi-protein FASTA from provided features"""
    all_fastas = []
    for feature in features:
        gene_id = get_gene_id(feature)
        fasta_seq = feature.qualifiers['translation'][0]

        # Never write empty fasta entries
        if len(fasta_seq) == 0:
            logging.debug("No translation for %s, skipping" % gene_id)
            continue

        all_fastas.append(">%s\n%s" % (gene_id, fasta_seq))
    full_fasta = "\n".join(all_fastas)
    return full_fasta
