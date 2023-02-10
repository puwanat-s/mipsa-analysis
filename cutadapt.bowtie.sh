gunzip -k *.fastq.gz

for file in $(ls | grep ".fastq" | sed 's/.fastq//'); do cutadapt --cores=8 -u -121 -o ${file}\_trimmed.fastq ${file}.fastq; done

for file in $(ls | grep ".fastq" | sed 's/_trimmed.fastq//'); do bowtie -p 8 -S -a -v 1 --best --strata -x ~/data-hlarman1/MIPSA-pep_db/RawData/MIPSApep_001/UCI.Pepseq.dictionary/UCI.PepSeq.dictionary.12142022 ${file}\_trimmed.fastq ${file}\.sam &> ${file}\.log; done

for i in *.log; do echo $i; echo; cat $i; echo; done > alignment_results.txt

rm *.log

r-studio-server.sh -n 1 -c 8 -t 2:0:0 -p defq

