# In contrast to the "traditional" yeast NET-seq protocol (Churchman 2012 - PMID 21248844), we used Bioo Scientific Small RNA-seq kit v3:
# - Reads were sequenced on PE mode;
# - Each read starts with a 4 nt UMI;
# - The first non-UMI base of R2 corresponds to the position of transcriptionally engaged RNAPII;
# For more details on our custom lib prep protocol, see Kindgren et al., 2019 (PMID 31863587);


# Download FASTQ files from our repository (observe that the data are PE):
# (wt, hda1-3)
# ..............
# ..............
# ..............


# Count raw reads:
for file in s*R2.fastq.gz; do 
  echo $file \
    $(( $(zcat $file | wc -l | awk '{print $1}') / 4 )); 
done

# Process UMIs (PE mode):
for f1 in s*R1.fastq.gz; do 
  f2=${f1/_R1/_R2} && 
  echo $f1 $f2 && 
  umi_tools extract --stdin=${f1} --read2-in=${f2} \
    --bc-pattern=NNNN --bc-pattern2=NNNN \
    --stdout=${f1/.fastq.gz/_UMI.fq.gz} \
    --read2-out=${f2/.fastq.gz/_UMI.fq.gz}; 
done

# Align to SacCer3 (SE mode, R2 only; observe the trimming of Illumina Small RNA-seq R2 adapter):
for f2 in s*R2_UMI.fq.gz; do 
  echo $f2 && 
  STAR --genomeDir saccer3_star --readFilesIn $f2 \
    --runThreadN 4 --outFileNamePrefix ${f2/R2_UMI.fq.gz/} \
    --outSAMmultNmax 1 --alignEndsType Extend5pOfRead1 \
    --readFilesCommand zcat \
    --clip3pAdapterSeq GATCGTCGGACTGTAGAACTCTGAAC \
    --outSAMtype BAM Unsorted; 
done

for file in *Aligned*; do mv $file ${file/_Aligned.out/}; done
rmdir *STARtmp; rm *out *tab

# Sort BAM files:
for file in s*bam; do 
  echo $file && 
  samtools sort $file -o ${file/.bam/_sorted.bam} && 
  rm $file; 
done

# Count aligned reads:
for file in s*bam; do 
  echo $file $(samtools flagstat $file | 
    sed -n '1p' | awk '{print $1}'); 
done 

# Remove mitochondrial, rRNA and tRNA reads:
for file in s*sorted.bam; do 
  echo $file && 
  bedtools intersect -v -abam $file \
    -b SacCer3_rRNA_tRNA_genes_ext100bp.bed > \
    ${file/.bam/_clean.bam}; 
done

# Remove multimappers:
for file in s*clean.bam; do 
  echo $file && 
  samtools view -hb -q 10 $file > ${file/.bam/_mapq.bam}; 
done

# Deduplicate (UMI-Tools):
for file in s*mapq.bam; do 
  echo $file && 
  samtools index $file && 
  umi_tools dedup --stdin=${file} \
    --stdout=${file/.bam/_dedup.bam}; 
done

# Make stranded Bedgraph files (with strand switch):
for str in "+" "-"; do 
  [ "$str" = "+" ] && n="rev" || n="fw"; 
  for file in s*dedup.bam; do 
    sample=${file/_sorted_clean_mapq_dedup.bam/} && 
    echo $n $sample && 
    bedtools genomecov -ibam $file -bg -5 \
      -strand $str > ${sample}_${n}.bg; 
  done; 
done

# Merge forward and reverse Bedgraph files: 
for f1 in s*fw.bg; do 
  f2=${file1/_fw/_rev} && 
  outf=${f1/_fw.bg/.bedgraph.gz} && 
  echo $f1 "+" $f2 "=" $outf && 
  awk 'BEGIN{OFS="\t"}{print $1,$2,$3,"-"$4}' $f2 | 
  cat $f1 - | sort -k1,1 -k2,2n | 
  sed '1i track type=bedGraph color=0,100,200 altColor=200,100,0' | 
  gzip > $outf; 
done

# Normalize to 1M tags:
for file in s*bedgraph.gz; do 
  norm=$( zcat $file | sed '/^[#t]/d' | \
    awk 'BEGIN{SUM=0}{SUM+=sqrt($4^2)*($3-$2)}\
      END{print SUM / 1000000}' ) && 
  echo $file $norm && 
  zcat $file | 
  awk -v norm=$norm 'BEGIN{OFS="\t"}{if ($0~/^[#t]/) \
    print $0; else print $1, $2, $3, $4 / norm}' | 
  gzip > ${file/.bedgraph.gz/_norm1M.bedgraph.gz}; 
done

