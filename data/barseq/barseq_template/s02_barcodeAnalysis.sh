#!/bin/bash                                                                                                                                                               

#SBATCH --job-name=barcodeAnalysis
#SBATCH --array=1
#SBATCH --mem=10G
#SBATCH --cpus-per-task=8
#SBATCH --output=slurm_out/barcode/z_bc_%a_%A.out
#SBATCH --error=slurm_out/barcode/z_bc_%a_%A.out
#SBATCH --mail-type=END
#SBATCH --mail-user=kevin.blake@wustl.edu

eval $( spack load --sh /duhk6u2 )
eval $( spack load --sh usearch )
eval $( spack load --sh py-umi-tools )

#A is forward read, B is reverse
#C is association of outer barcodes with sample
#D is association of inner barcodes with promoter
#E is output file
#F is list of identified barcodes
#G is input suffix to UMI-tools
basedir="$PWD"
indir="${basedir}/d01_clean_reads"
outidr="${basedir}/d02_barcodeanalysis"
sampledir="${basedir}/samplelists"

# read in sample name                                                                                                                                                     
sample=`sed -n ${SLURM_ARRAY_TASK_ID}p ${basedir}/samplelists/samplelist.txt`

A="${indir}/${sample}_FW_clean.fastq"
B="${indir}/${sample}_RV_clean.fastq"
C="${sampledir}/samplebarcodes.txt"
D="${sampledir}/genebarcodes.txt"
E="expression.txt"
F="barcodes.txt"
G="tabUMIs.txt"


set -x
time usearch -fastq_mergepairs $A -reverse $B \
        -threads 8 \
        -fastq_maxmergelen 94 -fastq_minmergelen 70 \
        -fastqout aligned.fastq

python barseq.py aligned.fastq $C $D $E $F $G

for file in *tabUMIs.txt
do
        umi_tools count_tab -I ${file} -S ${file}counts.txt
done

RC=$?
set +x
