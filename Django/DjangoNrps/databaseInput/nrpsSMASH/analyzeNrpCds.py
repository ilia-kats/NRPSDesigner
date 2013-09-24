import utils
import time
import random
from argparse import Namespace
from os import path
import os
import shutil
from hmmscanparser import parse_hmmscan_results
from substrates import run_nrps_substr_spec_predictions, calculate_consensus_prediction, generate_domainnamesdict, parse_nrps_preds

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Alphabet import IUPAC

from celeryHelper.helpers import update_celery_task_state_log

# main input
#withinclustergenes example
#[SeqFeature(FeatureLocation(ExactPosition(0), ExactPosition(3267), strand=1), type='CDS')]

#  decide later on regarding how to best set those options for server etc
#  possibly convert options from argparse.Namespace to dictionary
#  but for now use this as a hack.. WHICH has terrible side-effects


def nrpsSmash(dnaSeq):
    options = Namespace()
    options.outputfoldername = "/tmp/nrpspks_predictions_txt"
    options.record_idx = "" # used in NRPSPredictor2.nrpscodepred, check later what to set it to
    options.eukaryotic = 0
    tstFeature = SeqFeature(FeatureLocation(0, len(dnaSeq)), type="CDS", strand=1)
    tstFeature.qualifiers = {'gene':['gene']}
    sequenceRecord = SeqRecord(Seq(dnaSeq, IUPAC.unambiguous_dna),
        id = "sauce",
        name = "bbqSauce",
        description = "wtfDude")
    sequenceRecord.features = [tstFeature]
    analysis = specific_analysis(sequenceRecord, options)
    shutil.rmtree(options.raw_predictions_outputfolder)
    return analysis

def create_pksnrpsvars_object():
    pksnrpsvars = utils.Storage()
    pksnrpsvars.nrpspkstypedict = {}
    pksnrpsvars.compound_pred_dict = {}
    pksnrpsvars.domaindict = {}
    pksnrpsvars.failedstructures = []
    pksnrpsvars.dockingdomainanalysis = []
    pksnrpsvars.pksnrpscoregenes = []
    return pksnrpsvars


def create_rawoutput_storage_folder(options):
    if not os.path.exists(options.outputfoldername):
        os.mkdir(options.outputfoldername)
    #Make output folder for storing raw predictions
    options.raw_predictions_outputfolder = path.abspath(path.join(options.outputfoldername, str(time.time()) + str(random.random())))
    if not os.path.exists(options.raw_predictions_outputfolder):
        os.mkdir(options.raw_predictions_outputfolder)

def specific_analysis(seq_record, options):
    pksnrpsvars = create_pksnrpsvars_object()
    create_rawoutput_storage_folder(options)
    pksnrpsvars = annotate_pksnrps(pksnrpsvars, seq_record)
    if len(pksnrpsvars.pksnrpscoregenes) > 0:
        run_nrps_substr_spec_predictions(pksnrpsvars, seq_record, options)
        parse_nrps_preds(options, pksnrpsvars)
        calculate_consensus_prediction(pksnrpsvars)
        generate_domainnamesdict(pksnrpsvars)
    return pksnrpsvars

def annotate_pksnrps(pksnrpsvars, seq_record):
    withinclustergenes = seq_record.features #temporary hack...
    run_nrpspks_specific_hmmer(seq_record, withinclustergenes, pksnrpsvars)
    name_nrpspks(seq_record, pksnrpsvars, withinclustergenes)
    pksnrpsvars.pksnrpscoregenes = utils.get_pksnrps_cds_features(seq_record)
    return pksnrpsvars


def run_nrpspks_specific_hmmer(seq_record, withinclustergenes, pksnrpsvars):
	# lol this is a true multiFASTA protein file
    #nrpspksfasta = utils.get_specific_multifasta(seq_record, withinclustergenes)
    gene_id = "gene"
    fasta_seq = str(seq_record.seq.translate(to_stop=True))
    nrpspksfasta = ">%s\n%s" % (gene_id, fasta_seq)
    #antiSMASH actually checks for abMotifs here but remove for now, since Ive no idea what it does :P
   	#Analyse for C/A/PCP/E/KS/AT/ATd/DH/KR/ER/ACP/TE/TD/COM/Docking/MT/CAL domains
    # note from HMMER3 documentation: "TC thresholds are
	# generally considered to be the score of the lowest-scoring known true positive that
	# is above all known false positives."
    update_celery_task_state_log("-Scanning for NRP domains using HMMER3")
    nrpspksdomain_opts = ["--cut_tc"]
    nrpspksdomain_results = utils.run_hmmscan(utils.get_full_path(__file__, "nrpspksdomains.hmm"), nrpspksfasta, nrpspksdomain_opts)
    hmmlengthsdict = utils.hmmlengths(utils.get_full_path(__file__, "nrpspksdomains.hmm"))
    pksnrpsvars.domaindict = parse_hmmscan_results(nrpspksdomain_results, hmmlengthsdict)
    pksnrpsvars.domaindict2 = pksnrpsvars.domaindict
    filter_nonterminal_docking_domains(seq_record, pksnrpsvars)

def filter_nonterminal_docking_domains(seq_record, pksnrpsvars):
    dockingdomains = ['NRPS-COM_Nterm', 'NRPS-COM_Cterm']
    hitgenes = pksnrpsvars.domaindict.keys()
    feature_by_id = utils.get_feature_dict(seq_record)
    for hitgene in hitgenes:
        to_remove = []
        cdsfeature = feature_by_id[hitgene]
        cds_seq = str(seq_record.seq.translate(to_stop=True))
        hitgenelength = len(cds_seq)
        x = 0
        for hit in pksnrpsvars.domaindict[hitgene]:
            if hit[0] in dockingdomains:
                if not (hitgenelength - max(hit[1], hit[2]) < 50 or min(hit[1], hit[2]) < 50):
                    to_remove.append(x)
            x += 1
        to_remove.reverse()
        for idx in to_remove:
            del pksnrpsvars.domaindict[hitgene][idx]
        if pksnrpsvars.domaindict[hitgene] == []:
            del pksnrpsvars.domaindict[hitgene]


def name_nrpspks(seq_record, pksnrpsvars, withinclustergenes):
    pksnrpsvars.nrpspkstypedict = {}
    for feature in withinclustergenes:
        k = utils.get_gene_id(feature)
        if not pksnrpsvars.domaindict.has_key(k):
            continue
        if pksnrpsvars.domaindict[k] == []:
            continue
        #structure of domaindict: domaindict[genename] = [[name,start,end,evalue,score],[name,start,end,evalue,score], etc.]
        domainlist = []
        for i in pksnrpsvars.domaindict[k]:
            domainlist.append(i[0])

        if pksnrpsvars.domaindict.has_key(k):
            for i in pksnrpsvars.domaindict[k]:
                if ("Condensation_LCL" in domainlist or "Condensation_DCL" in domainlist or "Condensation_Starter" in domainlist or "Cglyc" in domainlist or "Condensation_Dual" in domainlist) and "AMP-binding" in domainlist and "PKS_KS" not in domainlist and "PKS_AT" not in domainlist:
                    nrpspkstype = "NRPS"
                elif ("Condensation_LCL" in domainlist or "Condensation_DCL" in domainlist or "Condensation_Starter" in domainlist or "Cglyc" in domainlist or "Condensation_Dual" in domainlist or "AMP-binding" in domainlist) and "PKS_KS" not in domainlist and "PKS_AT" not in domainlist:
                    nrpspkstype = "NRPS-like protein"
                else:
                    nrpspkstype = "the shit"
            if feature.qualifiers.has_key("sec_met"):
                feature.qualifiers['sec_met'].append("NRPS/PKS subtype: " + nrpspkstype)
            else:
                feature.qualifiers['sec_met'] = ["NRPS/PKS subtype: " + nrpspkstype]
            nrpspksdomains = pksnrpsvars.domaindict[k]

            for domain in nrpspksdomains:
                if feature.qualifiers.has_key("sec_met"):
                    feature.qualifiers['sec_met'].append("NRPS/PKS Domain: %s (%s-%s). E-value: %s. Score: %s;" % (domain[0], str(domain[1]), str(domain[2]), str(domain[3]), str(domain[4])))
                else:
                    feature.qualifiers['sec_met'] = ["NRPS/PKS Domain: %s (%s-%s). E-value: %s. Score: %s;" % (domain[0], str(domain[1]), str(domain[2]), str(domain[3]), str(domain[4]))]
        pksnrpsvars.nrpspkstypedict[k] = nrpspkstype
