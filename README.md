# handle RefSeq data
Get sequence metadata for RefSeq release catalogue using Entrez (also for suppressed or replaced sequences), filter for top 3 assemblies for each taxa and download Fastas 


- get_nuc_seqIds_taxIds.py: extract nucelotide sequence IDs from given RefSeq catalogue, and write sequence and taxonomy IDs out into 3 batches (-> seqId_taxId.batch{1,2,3}.txt.gz)

- download_seq_metadata.py: download metadata (seqId, seqLength, taxId, assemblyId) from Entrez for each sequence from given seqId_taxId file (use batches, distributed on different machines)

    $ python download_seq_metadata.py seqId_taxId.txt.gz seq_metadata.txt.gz 28 name@mail.com NCBI_API_key

- filter_ids.top3.py: load metadata, check, determine IDs of 3 longest assemblies for each taxa, and retrieve corresponding sequence IDs

- download_fastas.py: download fastas for selected ids
