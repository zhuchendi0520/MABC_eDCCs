ls *_1.fastq.gz|while read i
do a=${i%_1.fastq.gz}
b=${a}
j=${a}_2.fastq.gz
fq1=${a}_tr_1.fq
fq2=${a}_tr_2.fq
fq3=${a}_tr_S.fq
samp=${a}.paired.sam
sams=${a}.single.sam
bamp=${a}.paired.bam
bams=${a}.single.bam
bamm=${a}.merge.bam
sortbam=${a}.sort.bam
pileup=${a}.pileup
var=${a}.varscan
cns=${a}.cns
ppe=${a}.ppe
format=${a}.for
forup=${a}.forup
fix=${a}.fix
snp=${a}.snp
echo "sickle pe -l 35 -f $i -r $j -t sanger -o $fq1 -p $fq2 -s $fq3;bwa mem -t 1 -c 100 -R '@RG\\\tID:$b\\\tSM:$b\\\tPL:illumina' -M /home/zcd/script/NTM_mapping/MAB_mapping/ATCC19977.fasta $fq1 $fq2 > $samp;bwa mem -t 1 -c 100 -R '@RG\\\tID:$b\\\tSM:$b\\\tPL:illumina' -M /home/zcd/script/NTM_mapping/MAB_mapping/ATCC19977.fasta $fq3 > $sams;samtools view -bhSt /home/zcd/script/NTM_mapping/ATCC19977.fasta.fai $samp -o $bamp;samtools view -bhSt /home/zcd/script/NTM_mapping/ATCC19977.fasta.fai $sams -o $bams; samtools merge $bamm $bamp $bams;samtools sort $bamm -o $sortbam; samtools index $sortbam;samtools mpileup -q 30 -Q 20 -BOf /home/zcd/script/NTM_mapping/MAB_mapping/ATCC19977.fasta $sortbam > $pileup; java -jar /home/edwin/bin/VarScan.v2.3.6.jar mpileup2snp $pileup --min-coverage 3 --min-reads2 2 --min-avg-qual 20 --min-var-freq 0.01 --min-freq-for-hom 0.9 --p-value 99e-02 --strand-filter 0 > $var; java -jar /home/edwin/bin/VarScan.v2.3.6.jar mpileup2cns $pileup --min-coverage 3 --min-avg-qual 20 --min-var-freq 0.75 --min-reads2 2 --strand-filter 0 > $cns;perl /home/edwin/script/varscan_work_flow/0.1_PE_IS_filt_Rv.pl /home/zcd/script/NTM_mapping/M.abscessus_repeat_PEPPE_phage_transposase.ATCC19977.txt $var > $ppe;perl /home/edwin/script/varscan_work_flow/1_format_trans.pl $ppe > $format;perl /home/edwin/script/varscan_work_flow/2_fix_extract.pl $format > $fix;perl /home/edwin/script/varscan_work_flow/3.1_mix_pileup_merge.pl $format $pileup >$forup; cut -f2-4 $fix > $snp;rm $fq1 $fq2 $fq3 $samp $sams $bamp $bams $bamm $pileup $ppe $format $fix;"
done
