#!/bin/bash
#SBATCH --job-name="MIPSA_preprocessing_alignment"
#SBATCH --partition=parallel
#SBATCH -N 1
#SBATCH --cpus-per-task=8
#SBATCH --time=2:00:00
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=psangkh1@jhu.edu

ml cutadapt
ml bowtie
ml samtools
ml parallel

FASTQ_FILE=$(ls data/*.fastq.gz | head -n 1) # get one of the fastq files to calculate seq length
CPUS_PER_TASK=8

# Check if a fastq.gz file exists
if [[ -z "$FASTQ_FILE" ]]; then
    echo "No FASTQ files found in the directory!"
    exit 1
fi

echo "Using FASTQ file: $FASTQ_FILE"

# Extract the sequence length from the first read in the FASTQ file
SEQ_LEN=$(zcat "$FASTQ_FILE" | awk 'NR==2 {print length($0); exit}')

echo "Sequence length detected: $SEQ_LEN"

# Calculate the trimming value (-(SEQ_LEN - 61))
TRIM_VAL=$((-(SEQ_LEN - 61)))

# extract 41nt barcode from the reads
for file in $(ls data/*.fastq.gz | sed 's|data/||' | sed 's/_S[0-9]\+_R[12]_001\.fastq.gz//'); do cutadapt \
	--cores=$CPUS_PER_TASK \
	-u 20 \
	-u "$TRIM_VAL" \
	-o ${file}\_trimmed.fastq.gz data/${file}_S*_R*_001.fastq.gz; done

# align the barcodes with MIPSA dictionary (perfect match)
for file in $(ls data/*.fastq.gz | sed 's|data/||' | sed 's/_S[0-9]\+_R[12]_001\.fastq.gz//'); do bowtie \
	-p $CPUS_PER_TASK -S -a -v 0 \
	--sam-nosq \
	--best \
	--strata \
	-x dictionary/index_files ${file}\_trimmed.fastq.gz ${file}\.sam &> ${file}\.log; done


### create a tsv file containing stats from bowtie alignment
ALIGNMENT_OUTPUT="bowtie_results.tsv"

echo -e "samples\toriginal_reads\tmapped_reads\tpercent_alignment" > "$ALIGNMENT_OUTPUT" # add header

for log in *.log; do
    sample=$(basename "$log" .log)  # Extract sample name from filename

    # Extract relevant numbers from the log file
    original_reads=$(grep "# reads processed:" "$log" | awk '{print $4}')  # Total reads
    mapped_reads=$(grep "# reads with at least one alignment:" "$log" | awk '{print $8}')  # Mapped reads
    percent_alignment=$(grep "# reads with at least one alignment:" "$log" | awk -F '[()]' '{print $2}' | tr -d '%')  # Extracts percentage

    # Write data to TSV file
    echo -e "$sample\t$original_reads\t$mapped_reads\t$percent_alignment" >> "$ALIGNMENT_OUTPUT"

done

echo "Bowtie results saved in $ALIGNMENT_OUTPUT"

### convert sam to bam
for file in $(ls | grep ".sam" | sed 's/.sam//'); do samtools view \
	-@ $CPUS_PER_TASK -bt dictionary/MIPSA_UltimateHumanORF_dictionary.fa.fai \
	-o ${file}\.bam \
	${file}\.sam; done

rm *.sam
rm *.log
rm *_trimmed.fastq.gz

### generate count matrices from bam files
BAM_FILES=$(ls *.bam)

# Function to extract barcode counts
process_bam() {
    BAM_FILE=$1
    SAMPLE_NAME=$(basename "$BAM_FILE" .bam | awk -F'-' '{print $3"_"$4}') # choose your own indices for sample name and # of replicates
    OUTPUT_FILE="${SAMPLE_NAME}.tsv"

    # Extract read sequences, reference name, and filter only mapped reads (flag 0)
    samtools view "$BAM_FILE" | awk '$2 == 0 {print $10"-"$3}' | sort | uniq -c | awk '{print $2"\t"$1}' >> "$OUTPUT_FILE"

    echo "Processed: $BAM_FILE --> $OUTPUT_FILE"
}

export -f process_bam

parallel -j $CPUS_PER_TASK process_bam ::: $BAM_FILES

mkdir bowtie_output
mv *.bam bowtie_output
mv $ALIGNMENT_OUTPUT bowtie_output

# Sort all individual barcode count files by barcode
for f in *.tsv; do
    sort -k1,1 "$f" > "${f}.sorted"
done

cut -f1 *.tsv.sorted | sort -u > union.tsv

# Initialize the merged file with the union of barcodes
cp union.tsv merged_barcode_counts.tsv

# For each sorted file, join with the merged file
for f in *.tsv.sorted; do
    join -t $'\t' -a 1 -e 0 -o auto merged_barcode_counts.tsv "$f" > temp.tsv
    mv temp.tsv merged_barcode_counts.tsv
done

header="barcode"
for f in *.tsv.sorted; do
    sample=$(basename "$f" .tsv.sorted)
    header="${header}\t${sample}"
done
echo -e "$header" > header.tsv

# Prepend header to the merged file
cat header.tsv merged_barcode_counts.tsv > temp.tsv
mv temp.tsv merged_barcode_counts.tsv

# Clean up temporary files
rm header.tsv union.tsv *.tsv.sorted

echo "Barcode count matrix saved as merged_barcode_counts.tsv"

### create a study design file for enrichment analysis
STUDY_DESIGN="study_design.tsv"

echo -e "samples\tgroup" > "$STUDY_DESIGN" # header includes samples and group

# Extract full sample names (skipping the "barcode" column)
head -n1 merged_barcode_counts.tsv | cut -f2- | tr '\t' '\n' | \
awk -F'_' '{
    sample = $0;
    prefix = $1;
    if (!(prefix in grp)) {
        grp[prefix] = ++count;
    }
    print sample "\t" grp[prefix];
}' >> "$STUDY_DESIGN"

echo "Study design file created: $STUDY_DESIGN"

mkdir count_matrix
mv merged_barcode_counts.tsv study_design.tsv count_matrix
rm *.tsv