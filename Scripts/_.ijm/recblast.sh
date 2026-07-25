#!/bin/bash
#################################################################################################################
#SBATCH --job-name=sbatchTemplate ## Name of your job
#SBATCH --ntasks=16 ## number of cpu's to allocate for a job
#SBATCH --ntasks-per-node=16 ## number of cpu's to allocate per each node
#SBATCH --nodes=1 ## number of nodes to allocate for a job
#SBATCH --mem=128G ## memory to allocate for your job in MB
#SBATCH --time=1-00:00:00 ## time to allocate for your job in format: DD-HH:MM:SS
#SBATCH --error=%J.errors ## stderr file name(The %J will print job ID number)
#SBATCH --output=%J.output ## stdout file name(The %J will print job ID number)
#SBATCH --mail-type=NONE ## Send your job status via e-mail: Valid type values are NONE, BEGIN, END, FAIL, REQUEUE, ALL
########### Job information #############
echo "================================"
echo "Start at `date`"
echo "Job id is $SLURM_JOBID"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NTASKS processors."
echo "================================"
#########################################

######## Load required modules ##########
#. /etc/profile.d/modules.sh # Required line for modules environment to work
#module load openmpi/1.8.4 python/2.7 # Load modules that are required by your program
#conda init
#conda activate multiqc
#source /lustre1/home/mass/eskalon/miniconda/bin/activate agat
source /lustre1/home/mass/eskalon/miniconda/bin/activate blast
#########################################

### Below you can enter your program job command ###

#agat_convert_sp_gff2gtf.pl \
#-gff Stylophora_pistillata_gca002571385v1.Stylophora_pistillata_v1.60.gff3 \
#-o Stylophora_pistillata.ensembl.gtf

#makeblastdb \
  -in GCF_002571385.2_Stylophora_pistillata_v1.1_protein.faa \
  -dbtype prot \
  -parse_seqids \
  -out GCF_002571385.2_prot_db

#makeblastdb \
  -in biomin.genes.fasta \
  -dbtype prot \
  -parse_seqids \
  -out biomin_db

# forward search
blastp \
  -query biomin.genes.fasta \
  -db GCF_002571385.2_prot_db \
  -evalue 1e-5 \
  -max_target_seqs 20 \
  -seg yes \
 -num_threads 16 \
  -outfmt '6 qseqid sseqid pident length qlen slen qcovs qcovhsp evalue bitscore stitle' \
  -out biomin_vs_spis.tsv

# extract best 5 hits for each query
awk '$3 >= 25 && $7 >= 30 && $9 <= 1e-5' biomin_vs_spis.tsv |
sort -k1,1 -k10,10gr |
awk 'count[$1]++ < 5' \
> biomin_vs_spis.top5.filtered.tsv

# extract unique candidates
cut -f2 biomin_vs_spis.top5.filtered.tsv |
sort -u \
> candidate_spis_ids.txt


# extract target protein sequences
blastdbcmd \
  -db GCF_002571385.2_prot_db \
  -entry_batch candidate_spis_ids.txt \
  -out candidate_spis_hits.faa

# reverse search
blastp \
  -query candidate_spis_hits.faa \
  -db biomin_db \
  -evalue 1e-5 \
  -max_target_seqs 20 \
  -seg yes \
-num_threads 16 \
  -outfmt '6 qseqid sseqid pident length qlen slen qcovs evalue bitscore' \
  -out spis_vs_biomin.tsv

# best reverse hit
sort -k1,1 -k10,10gr spis_vs_biomin.tsv |
awk '!seen[$1]++' \
> spis_candidates_vs_biomin.best.tsv
