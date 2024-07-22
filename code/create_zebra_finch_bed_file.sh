#!/bin/bash
#### Extract the longest transcript

module load bioinfo-tools CGAT/0.3.3
source /sw/apps/bioinfo/CGAT/0.3.3/rackham/conda-install/etc/profile.d/conda.sh
conda activate base; conda activate cgat-s 

fasta_fai=data/external_raw/genome/Taeniopygia_guttata.taeGut3.2.4.dna.toplevel_final.fasta.fai

wget https://ftp.ensembl.org/pub/release-92/gtf/taeniopygia_guttata/Taeniopygia_guttata.taeGut3.2.4.92.gtf.gz

mkdir -p data/meta/

#### Selecting the longest transcript for each gene
zcat Taeniopygia_guttata.taeGut3.2.4.92.gtf.gz | cgat gtf2gtf --method=filter --filter-method=longest-transcript > data/meta/Taeniopygia_guttata.taeGut3.2.4.92.longestTranscript.gtf 

#### Get all CDS regions from each transcript (removing exons with 0 length)
cat data/meta/Taeniopygia_guttata.taeGut3.2.4.92.longestTranscript.gtf  | awk '$3=="exon" { print $0}' | awk '{print $1,$4,$5,$14}' | tr -d "\";" | sed 's/ /\t/g' | bedtools sort | awk '$2<$3 { print $0}' > data/meta/gene_coord_filenames_trans_exon.bed

#### Get all CDS regions from each transcript (removing exons with 0 length) with orientation
cat data/meta/Taeniopygia_guttata.taeGut3.2.4.92.longestTranscript.gtf | awk '$3=="exon" { print $0}' | awk '{print $1,$4,$5,$14,"dummy",$7}' | tr -d "\";" | sed 's/ /\t/g' | bedtools sort | awk '$2<$3 { print $0}' > data/meta/gene_coord_filenames_trans_exon_orientation.bed

#### Get all exon regions from each transcript (removing exons with 0 length), no overlap
cat data/meta/Taeniopygia_guttata.taeGut3.2.4.92.longestTranscript.gtf | awk '$3=="exon" { print $0}' | awk '{print $1,$4,$5,$14,"dummy",$7}' | tr -d "\";" | sed 's/ /\t/g' | bedtools sort | awk '$2<$3 { print $0}' | bedtools merge > data/meta/gene_coord_filenames_trans_exon_merge.bed

#### The whole gene region. This would contain the whole gene sequence and then we could always filter later
cat data/meta/Taeniopygia_guttata.taeGut3.2.4.92.longestTranscript.gtf  | awk '$3=="transcript" { print $0}' | awk '{print $1,$4,$5,$14}' | tr -d "\";" | sed 's/ /\t/g' | bedtools sort | awk '$2<$3 { print $0}' > data/meta/gene_coord_filenames_trans_whole_transcript.bed


#### Create a merged file so that freebayes does not call variants twice - use this one for freebayes
cat data/meta/Taeniopygia_guttata.taeGut3.2.4.92.longestTranscript.gtf  | awk '$3=="transcript" { print $0}' | awk '{print $1,$4,$5,$14}' | tr -d "\";" | sed 's/ /\t/g' | bedtools sort | awk '$2<$3 { print $0}' | bedtools merge > scratch/lost_W/gene_coord_filenames_trans_whole_transcript_merged.bed

#### Make a region file for freebayes-parallel
cat scratch/lost_W/gene_coord_filenames_trans_whole_transcript_merged.bed | awk '{print $1":"$2"-"$3}' > scratch/lost_W/gene_coord_filenames_trans_whole_transcript_merged_regions4freebayes.out

#### Make GFF file from GTF
gffread -E data/meta/Taeniopygia_guttata.taeGut3.2.4.92.longestTranscript.gtf > data/meta/Taeniopygia_guttata.taeGut3.2.4.92.longestTranscript.gff

### To get also intron as well as exonregions: 
cat $fasta_fai | cut -f 1-2 > data/external_raw/genome/Taeniopygia_guttata.taeGut3.2.4.dna.toplevel_final.genome.txt

bedtools sort -i data/meta/gene_coord_filenames_trans_exon.bed -g data/external_raw/genome/Taeniopygia_guttata.taeGut3.2.4.dna.toplevel_final.genome.txt | bedtools complement -i stdin -g data/external_raw/genome/Taeniopygia_guttata.taeGut3.2.4.dna.toplevel_final.genome.txt | bedtools intersect -a stdin -b data/meta/gene_coord_filenames_trans_whole_transcript.bed -wa -wb | cut -f 1-3,7 | awk -v OFS="\t" '{print $0,"intron"}' > data/meta/gene_coord_filenames_trans_exon_intron.bed

cat data/meta/gene_coord_filenames_trans_exon_intron.bed | bedtools sort -i stdin > data/meta/gene_coord_filenames_trans_exon_intron.sorted.bed

mkdir data/meta/exon_intron_ranges

# Lastly, exon and intron ranges used for the genome coverage calculations
cat data/meta/gene_coord_filenames_trans_exon_intron.sorted.bed | cut -f 4 | sort | uniq | while read trans ; do cat data/meta/gene_coord_filenames_trans_exon_intron.sorted.bed | grep $trans > data/meta/exon_intron_ranges/$trans.bed ; done

# And make gene list
cat data/meta/gene_coord_filenames_trans_exon_intron.bed | cut -f 4 | sort | uniq > data/meta/geneList.txt
### IMPORTANT ###

### New exon file with more info
zcat ../neosexchromosome/data/meta/Taeniopygia_guttata.taeGut3.2.4.92.gtf.gz | cgat gtf2gtf --method=filter --filter-method=longest-transcript | awk '$3=="CDS" { print }' | cut -f 1,4,5,7,9 | bedtools sort  | sed 's/;/\t/g' | cut -f 1-5,7,9 | tr -d "\";" | awk -v OFS="\t" '{print $6,$8,$10,$1,$2,$3}' > data/meta/zf.allChr.exon.list
