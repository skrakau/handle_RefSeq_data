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

seqIds_file = sys.argv[1]
file_out = sys.argv[2]
threads = int(sys.argv[3])
my_email = sys.argv[4]
my_api_key = sys.argv[5]
print("\nProcessing file: ", seqIds_file, flush=True)
print("Using ", threads, " threads.\n", flush=True)

with open(seqIds_file, "rt") as infile:
        selected_seqIds = [ row[0] for row in csv.reader(infile, delimiter = "\t") ]

##############################################################
# 4) retrieve fasta sequences for filtered sequence accessions
##############################################################
# Create a new downloader instance
edl2 = EntrezDownloader(
        email = my_email,
        api_key = my_api_key, 
        num_threads = threads,                   # avoid swp
        batch_size = 50,                        
        pbar = True                             
        )

import threading
tempfile_lock = threading.Lock()
def fasta_parser(data):
        ids = []
        with tempfile_lock:
                for rec in SeqIO.parse(io.StringIO(data), 'fasta'):
                        SeqIO.write(rec, handle_out, "fasta")  
                        ids.append(rec.id)   
                handle_out.flush()        
        return ids


print("Going to download ", len(selected_seqIds), "selected fasta sequences!", flush=True)
print("Writing results to file: ", file_out, flush=True)
handle_out =  open( file_out, "wt" )    # do not zip, make sure usable output in case of problem...
# Execute parallel efetch for the specified fastas
results, failed = edl2.efetch(
    db = 'nuccore',
    retmode = 'txt',
    rettype = 'fasta',
    ids = selected_seqIds,
    result_func = fasta_parser
)

if len(failed) > 0:
        print("    Download failed for ", len(failed), " sequences.", flush=True)
        pprint(failed)
        print("    Trying again with batchsize 1 ...", flush=True)
        edl2_failed = EntrezDownloader(
                email = my_email,
                api_key = my_api_key,
                num_threads = threads,
                batch_size = 10,
                pbar = True
        )
        temp_ids = failed
        for _ in range(100):
                temp_results, failed = edl2_failed.efetch(db = 'nuccore',ids = temp_ids,result_func = fasta_parser)
                results += temp_results
                if not failed:
                        break
                temp_ids = failed
        if failed:
                print("    still have ", len(failed), "failed results after multiple attempts.", flush=True)



handle_out.close()
print("Finished.", flush=True)

