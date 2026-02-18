configfile: "config.yaml"
SEED_FASTA=config['seed_fasta'] # "/home/daffa/Work/2025/11-JSPP67/hmm/seed/NP_390565.2.faa"
PREFIX=config['prefix'] # "prpF_homologs_in_plant"
RESULTS_DIR=config['results_dir'] # "/home/daffa/Work/2025/11-JSPP67/hmm/hmm_pipeline/prpF"

TARGET="/home/daffa/Work/2025/11-JSPP67/hmm/target/GCF_000003195.3_Sorghum_bicolor_NCBIv3_protein.faa"
LOCAL_BLASTDB_DIR="/home/daffa/Work/2025/11-JSPP67/hmm/hdd_link/blastdb/"
LOCAL_BLASTDB="/home/daffa/Work/2025/11-JSPP67/hmm/hdd_link/blastdb/refseq_protein"
TAXID_FILE="/home/daffa/Work/2025/11-JSPP67/taxdump/viridiplantae_taxids.txt"

BLAST_HITS=f"{RESULTS_DIR}/{PREFIX}-blast_hits.txt"
FILTERED_ACCESSIONS=f"{RESULTS_DIR}/{PREFIX}-accessions.txt"
FILTERED_SEQUENCES=f"{RESULTS_DIR}/{PREFIX}.faa"
ALIGNED_FILTERED_SEQUENCES=f"{RESULTS_DIR}/{PREFIX}.aln.faa"
HMM_DB=f"{RESULTS_DIR}/{PREFIX}.hmm"
PRESSED_HMM_DB=expand("{db}.{suffix}", db=HMM_DB, suffix=["h3f", "h3i", "h3m", "h3p"])
FINAL_RESULT=f"{RESULTS_DIR}/{PREFIX}-results.tbl"

'''
A sensitive homology detection strategy:
1. BLAST find close homologs.
2. HMM captures evolutionary patterns from these homologs. Finds genes in sorghum that have conserved amino acids regions as in the homologs.
'''

rule all:
    input: FINAL_RESULT
    # input: FILTERED_BLAST_HITS

rule blast_seed_gene:
    input: SEED_FASTA
    params:
        blastdb_dir = LOCAL_BLASTDB_DIR,
        db = LOCAL_BLASTDB,
        outfmt = "6 qacc sacc pident length evalue bitscore staxids",
        max_target_seqs = 10,
        taxids = 33090
    output: BLAST_HITS
    shell:
        '''
        export BLASTDB="{params.blastdb_dir}"
        blastp -db {params.db} -query {input} -outfmt "{params.outfmt}" -taxids {params.taxids} -max_target_seqs {params.max_target_seqs} > {output}
        '''

rule obtain_accessions:
    input: BLAST_HITS
    output: FILTERED_ACCESSIONS
    shell:
        '''
        cut -f2 {input} > {output}
        '''

rule download_sequences:
    input: FILTERED_ACCESSIONS
    output: FILTERED_SEQUENCES
    params:
        db = "protein"
    shell: 
        '''
        > {output}
        while read acc; do
            if [ ! -z "$acc" ]; then
                efetch -db {params.db} -id "$acc" -format fasta >> {output}
                sleep 0.1
            fi
        done < {input}
        '''

rule align_filtered_sequences:
    input: FILTERED_SEQUENCES
    output: ALIGNED_FILTERED_SEQUENCES
    shell: 
        '''
        mafft --auto {input} > {output} 
        '''

rule build_hmm_db:
    input: ALIGNED_FILTERED_SEQUENCES
    output: HMM_DB
    shell: 
        '''
        hmmbuild {output} {input}
        '''

rule hmmpress:
    input: HMM_DB
    output: PRESSED_HMM_DB
    shell:
        '''
        hmmpress {input}
        '''

rule hmmsearch:
    input:
        db = HMM_DB,
        target = TARGET,
        pressed = PRESSED_HMM_DB
    output: FINAL_RESULT
    shell:
        '''
        hmmsearch --tblout {output} {input.db} {input.target}
        '''
