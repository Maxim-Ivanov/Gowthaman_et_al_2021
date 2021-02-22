# There are 4 yeast NET-seq papers with WT samples;
# All of them used the same lib prep protocol (originally described in Churchman a. Weissman, 2011);
#
# 1) Churchman a. Weissman, 2011 (PMID 21248844)
# "Nascent transcript sequencing visualizes transcription at nucleotide resolution"
# 
# 2) Marquardt et al., 2014 (PMID 24949978)
# "A chromatin-based mechanism for limiting divergent noncoding transcription"
# 
# 3) Harlen et al., 2016 (PMID 27239037)
# "Comprehensive RNA Polymerase II Interactomes Reveal Distinct and Varied Roles for Each Phospho-CTD Residue"
#
# 4) Fischl et al., 2017 (PMID 28190769)
# "Paf1 Has Distinct Roles in Transcription Elongation and Differential Transcript Fate"


# Download FASTQ files from SRA (all data are SE):
echo -e "SRR072814\tNETseq_Churchman2012_wt_part1
SRR072815\tNETseq_Churchman2012_wt_part2
SRR072816\tNETseq_Churchman2012_wt_part3
SRR072817\tNETseq_Churchman2012_wt_part4
SRR072818\tNETseq_Churchman2012_wt_part5
SRR1197973\tNETseq_Marquardt2014_wt_rep1
SRR1197974\tNETseq_Marquardt2014_wt_rep2
SRR1197979\tNETseq_Marquardt2014_cac2_rep1
SRR1197980\tNETseq_Marquardt2014_cac2_rep2
SRR2005997\tNETseq_Harlen2016_wt_rep1
SRR2005998\tNETseq_Harlen2016_wt_rep2" > SacCer_NET-seq_acc.txt

awk '{cmd="fastq-dump --gzip "$1" && \
  mv "$1".fastq.gz "$2".fastq.gz"; \
  system(cmd)}' SacCer_NET-seq_acc.txt

ftp_prefix="ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR168"
out_prefix="NETseq_Fischl2017_wt"
wget ${ftp_prefix}/009/ERR1680009/ERR1680009.fastq.gz \
  -O ${out_prefix}_rep1.fastq.gz
wget ${ftp_prefix}/000/ERR1680010/ERR1680010.fastq.gz \
  -O ${out_prefix}_rep2.fastq.gz

# Concatenate FASTQ files corresponding to the same biological replicate in Churchman 2012:
zcat NETseq_Churchman2012_wt_part?.fastq.gz | 
gzip > NETseq_Churchman2012_wt.fastq.gz && 
rm NETseq_Churchman2012_wt_part?.fastq.gz

# Download SacCer3 genome from UCSC:
wget http://hgdownload.soe.ucsc.edu/goldenPath/sacCer3/bigZips/sacCer3.fa.gz
gunzip sacCer3.fa.gz

# Generate STAR index:
STAR --runMode genomeGenerate \
  --genomeFastaFiles sacCer3.fa \
  --runThreadN 4 --genomeDir saccer3_star

# Align to the SacCer3 genome (observe trimming the custom 3' adapter):
for file in NETseq*fastq.gz; do 
  echo $file && 
  STAR --genomeDir saccer3_star --readFilesIn $file \
    --runThreadN 4 --outFileNamePrefix ${file/.fastq.gz/_} \
    --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 \
    --readFilesCommand zcat --clip3pAdapterSeq ATCTCGTATGCCG \
    --outSAMtype BAM Unsorted; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rmdir *STARtmp; rm *out *tab

# Sort BAM files:
for file in NETseq*bam; do 
  echo $file && 
  samtools sort $file -o ${file/.bam/_sorted.bam}; 
done

# Download gene annotation for sacCer3:
wget http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/archive/saccharomyces_cerevisiae.20170114.gff.gz

# Create BED file with coordinates of known rRNA and tRNA genes:
zcat saccharomyces_cerevisiae.201701014.gff.gz | 
sed -n '1,/##FASTA/p' | 
sed '/^#/d;s/^chrmt/chrM/' | 
awk 'BEGIN{OFS="\t"}{if ($3=="rRNA_gene" || $3=="tRNA_gene") \
  print $1,$4-100,$5+100,$3,100,$7}' > \
  SacCer3_rRNA_tRNA_ext100bp.bed

# Add the whole chrM to the BED file:
echo -e "chrM\t1\t85779\tchrM_fw\t100\t+\nchrM\t1\t85779\tchrM_rev\t100\t-" >> SacCer3_rRNA_tRNA_ext100bp.bed

# Remove mitochondrial, rRNA and tRNA reads:
for file in NETseq*sorted.bam; do 
  echo $file && 
  bedtools intersect -v -abam $file \
    -b SacCer3_rRNA_tRNA_ext100bp.bed > \
    ${file/.bam/_clean.bam}; 
done

# Remove multimappers:
for file in NETseq*clean.bam; do 
  echo $file && 
  samtools view -hb -q 10 $file > ${file/.bam/_mapq.bam}; 
done

# Merge all wt BAM files into one super-sample:
samtools merge NETseq_merged_wt_sorted_clean_mapq.bam NETseq*_wt_*mapq.bam

# Optionally, remove the individual replicates:
# (except for wt and cac2 samples in Marquardt2014)
rm NETseq_Churchman2012*bam NETseq_Harlen2016*bam NETseq_Fischl2017*bam

# Make stranded Bedgraph files from either 5' terminal bases or whole NET-seq reads (with strand switch):
for str in "+" "-"; do 
  [ "$str" = "+" ] && 
  n="rev" || n="fw"; 
  for file in NETseq*mapq.bam; do 
    sample=${file/_sorted_clean_mapq.bam/} && 
    echo $n $sample && 
    bedtools genomecov -ibam $file -bg -5 \
      -strand $str > ${sample}_first_base_${n}.bg &&
    bedtools genomecov -ibam $file -bg \
      -strand $str > ${sample}_whole_read_${n}.bg; 
  done; 
done

# Merge forward and reverse Bedgraph files corresponding to the same sample:
for f1 in NETseq*fw.bg; do 
  f2=${f1/_fw/_rev} && 
  outf=${f1/_fw.bg/.bedgraph.gz} && 
  echo $f1 "+" $f2 "=" $outf && 
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $2 | 
  cat $f1 - | sort -k1,1 -k2,2n | 
  sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | 
  gzip > $outf; 
done

# Normalize "first_base" Bedgraph files to 1M tags:
for file in NETseq*first_base.bedgraph.gz; do 
  norm=$( zcat $file | sed '/^[#t]/d' | \
    awk 'BEGIN{SUM=0}{SUM+=sqrt($4^2)*($3-$2)}\
      END{print SUM / 1000000}' ) && 
  echo $file $norm && zcat $file | 
  awk -v norm=$norm 'BEGIN{OFS="\t"}{if ($0~/^[#t]/) print $0; \ 
    else print $1, $2, $3, $4 / norm}' | 
  gzip > ${file/.bedgraph.gz/_norm1M.bedgraph.gz}; 
done

# Use the NETseq_merged_wt_first_base_norm1M.bedgraph.gz file for visualization in the IGV browser (Fig. 1A, Fig. 3B-C, SFig. 1A-B, SFig. 3D-E);
# Also use this file to calculate FPKM values (required to stratify DNC transcripts by expression level on SFig. 3A-C);

# Use the NETseq_merged_wt_whole_read.bedgraph.gz for calling nascent-only transcripts (see 03-Correct_and_expand_SacCer3_annotation.R);

# Use the NETseq_Marquardt2014*rep[12]_first_base.bedgraph.gz files to calculate coding/non-coding ratios for all yeast genes (SFig. 1C).