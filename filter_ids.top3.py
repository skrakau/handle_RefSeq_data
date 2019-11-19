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
import re

from Bio import Entrez, SeqIO
from pprint import pprint
from datetime import datetime
from EntrezDownloader import EntrezDownloader
from collections import Counter

import sys

file_seqId_taxId = sys.argv[1]
file_ids = sys.argv[2]
print("\nProcessing file: ", file_ids, flush=True)


with gzip.open(file_seqId_taxId, "rt") as infile1:
        dict_seqId_taxId = { str(row[0]):str(row[1]) for row in csv.reader(infile1, delimiter = "\t") if row[0][:2] in ["NC", "NG", "NT", "NW", "NZ"] } 

with gzip.open(file_ids, "rt") as infile2:
        results = [ row for row in csv.reader(infile2, delimiter = "\t") ]


# test if entry is NCBI master record: NX_XXXX00000000.1 (representing a WGS project)
def is_master(seqId):
    return bool(re.search('[A-Z]{4}0+(\.\d){0,}$', seqId))
    

print("Read in ", len(results), " rows.", flush=True)
print("\nSet up dictionaries ...", flush=True)
# store taxId - assemblies dict
dict_taxId_assemblyIds = {}
dict_assemblyId_seqIds = {}
dict_seqId_seqLength = {}
for seqId, seqLength, taxId, assemblyId in results:
        # if not master record
        if not is_master(seqId):
                # sanity check
                if taxId != dict_seqId_taxId[seqId]:
                        print("Warning: sequence ", seqId, " had tax. ID ", dict_seqId_taxId[seqId], " assigned in RefSeq release catalogue, but got tax. ID ", taxId, " assigned in Entrez", flush=True)
                #assert taxId == dict_seqId_taxId[seqId],"taxonomy ID {taxId} wrong!"

                # update dict_taxId_assemblyIds
                if assemblyId == 'None':        # read in as str
                        assemblyId = "dummy_assembly_accession_taxid" + str(taxId)
                if taxId not in dict_taxId_assemblyIds:
                        dict_taxId_assemblyIds[taxId] = [assemblyId]
                else:
                        if assemblyId not in dict_taxId_assemblyIds[taxId]:
                                dict_taxId_assemblyIds[taxId].append(assemblyId)

                # update dict_assemblyId_seqIds
                if assemblyId not in dict_assemblyId_seqIds:
                        dict_assemblyId_seqIds[assemblyId] = [seqId]
                else:
                        if seqId not in dict_assemblyId_seqIds[assemblyId]:
                                dict_assemblyId_seqIds[assemblyId].append(seqId)                            

                # update dict_seqId_seqLength
                dict_seqId_seqLength[seqId] = int(seqLength)

######################################################
# 2) get assembly lengths and keep only 3 largest ones
######################################################
print("Compute assembly lengths and select top 3 assemblies for each taxa...", flush=True)


# compute assembly lengths from sequence lengths
dict_assemblyId_assemblyLength = {}
for assemblyId in dict_assemblyId_seqIds.keys():
        length = 0
        if assemblyId[:5] != "dummy":
                for seqId in dict_assemblyId_seqIds[assemblyId]:
                        length += dict_seqId_seqLength[seqId]

        dict_assemblyId_assemblyLength[assemblyId] = length


# for each taxId select top 3 assemblies
# (sequences assigned to no assembly will only be used, if less than 3 assemblies exist for taxId)
selected_seqIds = []
for taxId in dict_taxId_assemblyIds.keys():
        # get subset of taxId specific assemblyIds - assemblyLength
        current_assemblyId_lengths = {i: dict_assemblyId_assemblyLength.get(i, None) for i in dict_taxId_assemblyIds[taxId]}
        k = Counter(current_assemblyId_lengths)
        top3 = k.most_common(3) 
        for i in top3:
                selected_seqIds += dict_assemblyId_seqIds[i[0]]

#print([k for (k,v) in Counter(selected_seqIds).items() if v > 1])

# TODO check if length is same as reported on ncbi!!! (when using whole catalogue)
print("Check if assembly lengths are the same as reported on ncbi:", flush=True)
count = 0
for assemblyId in dict_assemblyId_assemblyLength.keys():
        print("    ", assemblyId, " :", dict_assemblyId_assemblyLength[assemblyId], " ", flush=True)
        count += 1
        if count == 10:
                break


########################################################################
# 3) write out filtered sequence accessions in case something goes wrong
########################################################################

print("For top3 assemblies keep ", len(selected_seqIds), " sequences.", flush=True)

id_file = open("selectd_seqIds.txt",'w')
for seqId in selected_seqIds:
        print(f"{seqId}", file = id_file, flush=True)

id_file.close()
