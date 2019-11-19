#!/usr/bin/env python3
####################################################################################################
#
# Author: Sabrina Krakau
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
####################################################################################################

import sys
import gzip
import csv
import xml.etree.ElementTree as ET
import tqdm
import io
import argparse

from Bio import Entrez, SeqIO
from pprint import pprint
from datetime import datetime
from EntrezDownloader import EntrezDownloader
from collections import Counter

import sys

#catalog = '/Users/skrakau/Development/RefSeq-release70.catalog.test.gz'         # test: 1000 sequences
file_in = sys.argv[1]   # seqId taxId
file_out = sys.argv[2]
threads = int(sys.argv[3])
my_email = sys.argv[4]
my_api_key = sys.argv[5]
print("\nProcessing file: ", file_in, flush=True)
print("Using ", threads, " threads.\n", flush=True)

with gzip.open(file_in, "rt") as infile:
        dict_seqId_taxId = { str(row[0]):str(row[1]) for row in csv.reader(infile, delimiter = "\t") if row[0][:2] in ["NC", "NG", "NT", "NW", "NZ"] }
        # NC_, NG_, NT_, NW_, NZ_ -> use only DNA sequences (should be already filtered)
        # add AC_ ?

seqIds = [ key for key in dict_seqId_taxId.keys() ]

##############################################
# 1) get assembly accessions for each sequence
##############################################
# Create a new downloader instance
edl = EntrezDownloader(
        email = my_email,                  # An email address. You might get blocked by the NCBI without specifying one.
        api_key = my_api_key,              # An API key. You can obtain one by creating an NCBI account. Speeds things up.
        num_threads = threads,             # The number of parallel requests to make
        batch_size = 50,                   # The number of IDs to fetch per request
        pbar = True                        # Enables a progress bar, requires tqdm package
        )


def get_tax_id(rec):
        try:
                f_table  = rec['GBSeq_feature-table']
                source   = [ elem for elem in f_table if elem['GBFeature_key']=='source' ][0]
                quals    = source['GBFeature_quals']
                db_xref  = [ elem for elem in quals if elem['GBQualifier_name']=='db_xref' and elem['GBQualifier_value'][:5]=='taxon'][0]
                tax_id   = int(db_xref['GBQualifier_value'].split(':')[1])
                return str(tax_id)
        except (ValueError, KeyError, IndexError):
                return None
    

def get_seq_id(rec):
        try:
                return str(rec['GBSeq_accession-version'])
        except KeyError:
                return None


def get_seq_length(rec):
        try:
                return int(rec['GBSeq_length'])
        except KeyError:
                return None


def get_assembly_id(rec):
        try:
                x_refs = rec['GBSeq_xrefs']
                assembly = [ elem for elem in x_refs if elem['GBXref_dbname']=='Assembly' ][0]
                return str(assembly['GBXref_id'])
        except (KeyError, IndexError, ValueError):
                return None
    

tmp_file = gzip.open(file_out,'wt')

import threading
tempfile_lock = threading.Lock()
def extract_ids(xml_text):
        batch_results =  []
        with tempfile_lock:
                for rec in Entrez.parse(io.StringIO(xml_text)):
                        seqId = get_seq_id(rec)
                        seqLen = get_seq_length(rec)
                        taxId = get_tax_id(rec)
                        assemblyId = get_assembly_id(rec)
                        print(f"{seqId}\t{seqLen}\t{taxId}\t{assemblyId}", file = tmp_file, flush=True)
                        batch_results.append((seqId, seqLen, taxId, assemblyId))
        return batch_results


# Execute parallel efetch for the specified IDs
print("Going to process ", len(seqIds), " RefSeq nucleotide sequences.", flush=True)
results, failed = edl.efetch(
        db = 'nuccore',
        ids = seqIds,
        #result_func = lambda xml_text : [ rec for rec in Entrez.parse(io.StringIO(xml_text)) ]
        result_func = extract_ids
)

if failed:
        print("    Access failed for ", len(failed), " sequences!", flush=True)
        pprint(failed)
        print("    Trying again with batchsize 1 ...", flush=True)
        edl_failed = EntrezDownloader(
                email = my_email,
                api_key = my_api_key,
                num_threads = threads,                 # avoid to many requests per second
                batch_size = 10,
                pbar = True
        )
        temp_ids = failed
        for _ in range(100):
                temp_results, failed = edl_failed.efetch(db = 'nuccore',ids = temp_ids,result_func = extract_ids)
                results += temp_results
                if not failed:
                        break
                temp_ids = failed
        if failed:
                print("    still have ", len(failed), "failed results after multiple attempts.", flush=True)


tmp_file.close()
print("Finished downloading of metainfo.", flush=True)

