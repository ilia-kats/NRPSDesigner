
# Copyright (C) 2010-2012 Marnix H. Medema
# University of Groningen
# Department of Microbial Physiology / Groningen Bioinformatics Centre
#
# Copyright (C) 2011,2012 Kai Blin
# University of Tuebingen
# Interfaculty Institute of Microbiology and Infection Medicine
# Div. of Microbiology/Biotechnology
#
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import logging
from substrates_nrps import *

import os
import utils

def run_nrps_substr_spec_predictions(pksnrpsvars, seq_record, options):
    #Predict NRPS A domain specificities with NRPSPredictor and Minowa et al. method

    # first extract A domains
    pksnrpsvars.nrpsnames, pksnrpsvars.nrpsseqs = extract_nrps_genes(pksnrpsvars.pksnrpscoregenes, pksnrpsvars.domaindict, seq_record, extra_aa=120)
    if len(pksnrpsvars.nrpsnames) > 0:
        run_nrpspredictor(seq_record, pksnrpsvars.nrpsnames, pksnrpsvars.nrpsseqs, options)
        run_minowa_predictor_nrps(pksnrpsvars.pksnrpscoregenes, pksnrpsvars.domaindict, seq_record, options)

def parse_nrps_preds(options, pksnrpsvars):
    pksnrpsvars.minowa_nrps_preds = {}
    pksnrpsvars.minowa_nrps_preds_details = {}
    pksnrpsvars.nrps_svm_preds = {}
    pksnrpsvars.nrps_svm_preds_details = {}
    pksnrpsvars.nrps_code_preds = {}
    pksnrpsvars.nrps_code_preds_details = {}
    substratetransdict2 = {'pipecolate':'pip','fOHOrn':'orn','beta-Lys':'blys','5NhOrn':'orn','OHOrn':'orn','Aad':'Aaa','bOHTyr':'bht'}
    if len(pksnrpsvars.nrpsnames) > 0:
        minowa_a_file = open(options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_minowa_nrpspredoutput.txt","r")
        minowa_a_file = minowa_a_file.read()
        minowa_a_file = minowa_a_file.replace("\r","\n")
        parts = minowa_a_file.split("\\\\\n")[1:]
        for i in parts:
            partlines = i.split("\n")
            acc = partlines[0]
            tophit = partlines[2].split("\t")[0]
            if tophit in substratetransdict2.keys():
                tophit = substratetransdict2[tophit]
            pksnrpsvars.minowa_nrps_preds[acc] = tophit.lower()
            pksnrpsvars.minowa_nrps_preds_details[acc] = "<b>Minowa HMM method A-domain<br>Substrate specificity prediction top hits:</b><br>\n" + partlines[1] + "<br>\n" + partlines[2] + "<br>\n" + partlines[3] + "<br>\n" + partlines[4] + "<br><br>\n\n"
        nrpspredictorfile1 = open(options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_nrpspredictor2_svm.txt","r")
        nrpspredictorfile2 = open(options.raw_predictions_outputfolder + os.sep + "ctg" + str(options.record_idx) + "_nrpspredictor2_codes.txt","r")
        nrpspredictorfile1 = nrpspredictorfile1.read()
        nrpspredictorfile1 = nrpspredictorfile1.replace("\r","\n")
        lines = nrpspredictorfile1.split("\n")[1:-1]
        for k in lines:
            tabs = k.split("\t")
            if tabs[6] != "N/A":
                pksnrpsvars.nrps_svm_preds[tabs[0]] = tabs[6]
            elif tabs[5] != "N/A":
                pksnrpsvars.nrps_svm_preds[tabs[0]] = tabs[5]
            elif tabs[4] != "N/A":
                pksnrpsvars.nrps_svm_preds[tabs[0]] = tabs[4]
            else:
                pksnrpsvars.nrps_svm_preds[tabs[0]] = tabs[3]
            pksnrpsvars.nrps_svm_preds_details[tabs[0]] = "<b> NRPSPredictor2 SVM prediction details:</b><br>\n8 Angstrom 34 AA code:<br>\n" + tabs[1] + "<br>\nPredicted physicochemical class:<br>\n" + tabs[3] + "<br>\nLarge clusters prediction:<br>\n" + tabs[4] + "<br>\nSmall clusters prediction:<br>\n" + tabs[5] + "<br>\nSingle AA prediction:<br>\n" + tabs[6] + "<br><br>\n\n"
        nrpspredictorfile2 = nrpspredictorfile2.read()
        nrpspredictorfile2 = nrpspredictorfile2.replace("\r","\n")
        lines = nrpspredictorfile2.split("\n")[:-1]
        for k in lines:
            tabs = k.split("\t")
            pksnrpsvars.nrps_code_preds[tabs[0]] = tabs[1]
            pksnrpsvars.nrps_code_preds_details[tabs[0]] = "<b> NRPSPredictor2 Stachelhaus code prediction:</b><br>\n" + tabs[1] + "<br><br>\n\n"

#not in db yet
# aaa
# mpro
# dhb
# pgly
# aeo
#4mha
#pico
#phg
# dha
# scy
# pip
#bmt
#adds
#HIV VS 2HIVA????
#bht
#4pPro
#ala-b
#dht
#sal
#tcl
#hyv-d
#iva
#vol
#mal
#mmal
#mxmal
#emal

def calculate_consensus_prediction(pksnrpsvars):
    #Combine substrate specificity predictions into consensus prediction
    pksnrpsvars.consensuspreds = {}
    available_smiles_parts = ['GLY','ALA','VAL','LEU','ILE','MET','PRO','PHE','TRP','SER','THR','ASN','GLN','TYR','CYS','LYS','ARG',
    'HIS','ASP','GLU','MPRO','ORN','PGLY','DAB','BALA','AEO','DHA','PIP','BMT','gly','ala','val','leu','ile','met','pro','phe','trp','ser',
    'thr','asn','gln','tyr','cys','lys','arg','his','asp','glu','aaa','mpro','dhb','2hiva','orn','pgly','dab','bala','aeo','4mha','pico','phg',
    'dha','scy','pip','bmt','adds','aad','abu','hiv','dhpg','bht','3-me-glu','4pPro','ala-b','ala-d','dht','Sal','tcl','lys-b','hpg','hyv-d',
    'iva','vol','mal','mmal','mxmal','emal','nrp','pk','Gly','Ala','Val','Leu','Ile','Met','Pro','Phe','Trp','Ser','Thr','Asn','Gln','Tyr',
    'Cys','Lys','Arg','His','Asp','Glu','Mpro','23Dhb','34Dhb','2Hiva','Orn','Pgly','Dab','Bala','Aeo','4Mha','Pico','Aaa','Dha','Scy','Pip',
    'Bmt','Adds','DHpg','DHB','nrp','pk']

    for feature in pksnrpsvars.pksnrpscoregenes:
        locus = utils.get_gene_id(feature)
        nra = 0
        j = pksnrpsvars.domaindict[locus]
        for k in j:
            if k[0] == "AMP-binding" or k[0] == "A-OX":
                nra +=1
                preds = []
                preds.append(pksnrpsvars.minowa_nrps_preds[locus + "_A" + str(nra)])
                preds.append(pksnrpsvars.nrps_svm_preds[locus + "_A" + str(nra)])
                preds.append(pksnrpsvars.nrps_code_preds[locus + "_A" + str(nra)])
                cpred = "n"
                for l in preds:
                    if preds.count(l) > 1:
                        if l in available_smiles_parts:
                            pksnrpsvars.consensuspreds[locus + "_A" + str(nra)] = l
                        else:
                            pksnrpsvars.consensuspreds[locus + "_A" + str(nra)] = "nrp"
                        cpred = "y"
                if cpred == "n":
                    pksnrpsvars.consensuspreds[locus + "_A" + str(nra)] = "nrp"


def generate_domainnamesdict(pksnrpsvars):
    pksnrpsvars.domainnamesdict = {}
    for feature in pksnrpsvars.pksnrpscoregenes:
        locus = utils.get_gene_id(feature)
        j = pksnrpsvars.domaindict[locus]
        domainnames = []
        for k in j:
            domainnames.append(k[0])
        pksnrpsvars.domainnamesdict[locus] = domainnames
