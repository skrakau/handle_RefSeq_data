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
import io
import csv
import argparse

#catalog = '/Users/skrakau/Development/RefSeq-release70.catalog.test.gz'         # test: 1000 sequences
catalog = sys.argv[1]
print("\nProcessing RefSeq catalog file: ", catalog, flush=True)
with gzip.open(catalog, "rt") as infile:
        dict_seqId_taxId = { str(row[2]):str(row[0]) for row in csv.reader(infile, delimiter = "\t") if row[2][:2] in ["NC", "NG", "NT", "NW", "NZ"] }
        # NC_, NG_, NT_, NW_, NZ_ -> use only DNA sequences

seqIds = [ key for key in dict_seqId_taxId.keys() ]
print("Going to process ", len(seqIds), " RefSeq nucleotide sequences.", flush=True)

# write out 3 batches
batch_size = int(len(seqIds)/3)

id_file1 = gzip.open("seqId_taxId.batch1.txt.gz",'wt')
for seqId in seqIds[:batch_size]:
        print(f"{seqId}\t{dict_seqId_taxId[seqId]}", file = id_file1, flush=True)
id_file1.close()

id_file2 = gzip.open("seqId_taxId.batch2.txt.gz",'wt')
for seqId in seqIds[batch_size:2*batch_size]:
        print(f"{seqId}\t{dict_seqId_taxId[seqId]}", file = id_file2, flush=True)
id_file2.close()

id_file3 = gzip.open("seqId_taxId.batch3.txt.gz",'wt')
for seqId in seqIds[2*batch_size:]:
        print(f"{seqId}\t{dict_seqId_taxId[seqId]}", file = id_file3, flush=True)
id_file3.close()

print("Finished filtering nucleotide sequences and creating batches.")