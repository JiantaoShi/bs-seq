import argparse, os

parser = argparse.ArgumentParser()
parser.add_argument('-I', "--input_bam", required=True,  help="Input BAM file.")
parser.add_argument('-S', "--snp_locus", required=True,  help="Genomic position of one SNP [chr:pos].")
parser.add_argument('-R', "--snp_ref", required=True,  help="SNP reference allele.")
parser.add_argument('-A', "--snp_alt", required=True,  help="SNP alternate allele.")
parser.add_argument('-Q', "--MAPQ", type=int, default=10, required=False, help="Minimum mapping quality of reads to be considered.")
parser.add_argument('-C', "--CpG", required=False, help="Indexed CpG annotation file for mHapSuite.")
parser.add_argument("-mHap", "--mHap", help="Convert resulting BAM to mHap files and calculate summary statistics.", action="store_true")
parser.add_argument("-debug", "--debug", default = False, help="Print bash command for debug.", action="store_true")
parser.add_argument('-O',  "--outFolder", required=True, help="Output folder.")
parser.add_argument('-tag',  "--tag", required=True, help="Output prefix.")
args = parser.parse_args()

def subprocess_wrap(cmd, debug):
    if debug:
        print(cmd)
        return
    os.system(cmd)

# bin
match_maker = '/sibcb1/bioinformatics/shijiantao/wgbsTk/wgbs_tools/src/pipeline_wgbs/match_maker'
snp_patter = '/sibcb1/bioinformatics/shijiantao/wgbsTk/wgbs_tools/src/pipeline_wgbs/snp_patter'
mHapSuite = '/sibcb1/bioinformatics/shijiantao/wgbsTk/mHapSuite-2.1-jar-with-dependencies.jar'

# output files
snp_chr, snp_pos = args.snp_locus.split(':')
out_ref = f'{args.outFolder}/{args.tag}.{snp_chr}.{snp_pos}.{args.snp_ref}'
out_alt = f'{args.outFolder}/{args.tag}.{snp_chr}.{snp_pos}.{args.snp_alt}'

# reads with reference allele
cmd = f'cat <(samtools view -H {args.input_bam}) <(samtools view -F 1797 -q {args.MAPQ} {args.input_bam} {snp_chr}:{snp_pos}-{snp_pos}'
cmd += f' | {match_maker} - | {snp_patter} --snp_pos {snp_pos} --snp_let1 {args.snp_ref} --snp_let2 {args.snp_alt})'
cmd += f' | samtools view -hb - | samtools sort -O bam - > {out_ref}.bam'
cmd += f' && samtools index {out_ref}.bam'
final_cmd = f'/bin/bash -c "{cmd}"'
subprocess_wrap(final_cmd, args.debug)

# reads with alternate allele
cmd = f'cat <(samtools view -H {args.input_bam}) <(samtools view -F 1797 -q {args.MAPQ} {args.input_bam} {snp_chr}:{snp_pos}-{snp_pos}'
cmd += f' | {match_maker} - | {snp_patter} --snp_pos {snp_pos} --snp_let1 {args.snp_alt} --snp_let2 {args.snp_ref})'
cmd += f' | samtools view -hb - | samtools sort -O bam - > {out_alt}.bam'
cmd += f' && samtools index {out_alt}.bam'
final_cmd = f'/bin/bash -c "{cmd}"'
subprocess_wrap(final_cmd, args.debug)


# summary
region = snp_chr + ':' + str(int(snp_pos) - 250) + '-' + str(int(snp_pos) + 250)
if args.mHap and args.CpG is not None:
	# convert BAM to mHap
	cmd = f'java -jar {mHapSuite} convert --cpgPath {args.CpG} --inputFile {out_ref}.bam --region {region} --outPutFile {out_ref}.mhap.gz'
	subprocess_wrap(cmd, args.debug)
	cmd = f'java -jar {mHapSuite} convert --cpgPath {args.CpG} --inputFile {out_alt}.bam --region {region} --outPutFile {out_alt}.mhap.gz'
	subprocess_wrap(cmd, args.debug)
	# index mHaps
	cmd = f'tabix -b 2 -e 3 {out_ref}.mhap.gz && tabix -b 2 -e 3 {out_alt}.mhap.gz'
	subprocess_wrap(cmd, args.debug)
	# region-level summary
	cmd = f'java -jar {mHapSuite} stat  --cpgPath {args.CpG} --mhapPath {out_ref}.mhap.gz'
	cmd += f' --metrics MM PDR CHALM MHL MCR MBS Entropy R2 --r2Cov 10 --region {region} --outputFile {out_ref}.tsv'	
	subprocess_wrap(cmd, args.debug)
	cmd = f'java -jar {mHapSuite} stat  --cpgPath {args.CpG} --mhapPath {out_alt}.mhap.gz'
	cmd += f' --metrics MM PDR CHALM MHL MCR MBS Entropy R2 --r2Cov 10 --region {region} --outputFile {out_alt}.tsv'
	subprocess_wrap(cmd, args.debug)
