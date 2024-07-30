import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-R1', "--Fastq_R1", required=True,  help="Read 1 fastq file.")
parser.add_argument('-R2', "--Fastq_R2", required=False, help="Read 2 fastq file.")
parser.add_argument('-N', "--threadN", type=int, default=8, required=False, help="Number of threads.")
parser.add_argument('-I', "--index", required=True,  help="BISCUIT index.")
parser.add_argument('-R', "--reference", required=True,  help="Reference in FASTA format.")
parser.add_argument("-rrbs", "--rrbs", help="whether it's RRBS data", action="store_true")
parser.add_argument("-runQC", "--runQC", help="whether run QC script", action="store_true")
parser.add_argument('-O',  "--outFolder", required=True, help="Output folder.")
parser.add_argument('-tag',  "--tag", required=True, help="Output prefix.")
args = parser.parse_args()

import sys,os,glob,re, uuid, subprocess

if not os.path.exists(args.outFolder):
    os.makedirs(args.outFolder)

biscuit = '/sibcb1/bioinformatics/shijiantao/miniforge3/bin/biscuit'

def subprocess_wrap(cmd, debug):
    if debug:
        print(cmd)
    os.system(cmd)

# alignment and deduplication
prefix_out = f'{args.outFolder}/{args.tag}'

cmd = f'{biscuit} align -@ {args.threadN} {args.index} {args.Fastq_R1}'
if args.Fastq_R2 is not None:
	if args.rrbs:
		cmd += f' {args.Fastq_R2}'
	else:
		cmd += f' {args.Fastq_R2} | /sibcb/program/bin/dupsifter {args.reference}'
cmd += f' | samtools view -b -o {prefix_out}_dup.bam'
subprocess_wrap(cmd, True)

# sort
bam_out = f'{prefix_out}.bam'
cmd = f'samtools sort -@ {args.threadN} -o {bam_out} {prefix_out}_dup.bam && samtools index {bam_out}'
subprocess_wrap(cmd, True)

# QC
if args.runQC:
	cmd = f'{biscuit} qc {args.reference} {bam_out} {prefix_out}'
	subprocess_wrap(cmd, True)

# Extract methylation
mCall = args.outFolder + '/mCall'
os.makedirs(mCall)
prefix_mCall = f'{mCall}/{args.tag}'
## pilepu
cmd = f'{biscuit} pileup -@ {args.threadN} -o {prefix_mCall}.vcf {args.reference} {bam_out}'
subprocess_wrap(cmd, True)
## index pilepu
cmd = f'bgzip -@ {args.threadN} {prefix_mCall}.vcf && tabix -p vcf {prefix_mCall}.vcf.gz'
subprocess_wrap(cmd, True)
## methylation call
cmd = f'{biscuit} vcf2bed -k 1 -t cg {prefix_mCall}.vcf.gz | {biscuit} mergecg -c {args.reference} -'
cmd += " | awk -v OFS='\t' '{ print $1, $2+1, $3-1, $4, $5, $6 }' >"
cmd += f' {prefix_mCall}_CpG.cov'
subprocess_wrap(cmd, True)

# clean
cmd = f'gzip {prefix_mCall}_CpG.cov && rm {prefix_out}_dup.bam'
subprocess_wrap(cmd, True)
