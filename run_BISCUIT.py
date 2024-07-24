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

import sys,os,glob,re, uuid

if not os.path.exists(args.outFolder):
    os.makedirs(args.outFolder)

biscuit = '/sibcb1/bioinformatics/shijiantao/miniforge3/bin/biscuit'
MethylDackel = '/sibcb3/bioinformatics3/shijiantao/miniforge3/bin/MethylDackel'

# alignment and deduplication
prefix_out = f'{args.outFolder}/{args.tag}'

cmd = f'{biscuit} align -@ {args.threadN} {args.index} {args.Fastq_R1}'
if args.Fastq_R2 is not None:
	if args.rrbs:
		cmd += f' {args.Fastq_R2}'
	else:
		cmd += f' {args.Fastq_R2} | dupsifter {args.reference}'
cmd += f' | samtools view -b -o {prefix_out}_dup.bam'
print(cmd)
os.system(cmd)

# sort
bam_out = f'{prefix_out}.bam'
cmd = f'samtools sort -@ {args.threadN} -o {bam_out} {prefix_out}_dup.bam'
print(cmd)
os.system(cmd)

# QC
if args.runQC:
	cmd = '/sibcb1/bioinformatics/shijiantao/miniforge3/bin/biscuit qc'
	cmd += f' {args.reference} {bam_out} {prefix_out}'
	print(cmd)
	os.system(cmd)

# Extract methylation
mCall = args.outFolder + '/mCall/'
if not os.path.exists(mCall):
    os.makedirs(mCall)
opref = mCall + args.tag
cmd = f'{MethylDackel} extract {args.reference} {bam_out} --mergeContext --opref {opref}'
print(cmd)
os.system(cmd)
cmd = f'{MethylDackel} extract {args.reference} {bam_out} --opref {opref} --CHG'
print(cmd)
os.system(cmd)
cmd = f'{MethylDackel} extract {args.reference} {bam_out} --opref {opref} --CHH'
print(cmd)
os.system(cmd)

# clean
cmd = f'rm {prefix_out}_dup.bam'
os.system(cmd)
