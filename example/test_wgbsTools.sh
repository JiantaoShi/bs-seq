python run_wgbsTools.py -I wgbs_tools/tutorial/bams/Left_Ventricle_STL001.IGF2.bam \
	-S chr11:2019496 -R C -A A \
	--CpG hg19_CpG.gz \
	--mHap -O . \
	--tag Left_Ventricle_STL001

python run_snpCov.py \
	-I /sibcb1/bioinformatics/shijiantao/WGBS_Normal/alignment/Z0000043W/Z0000043W.bam \
	-S out/dbsnp_146.hg38.common_nonCpG.vcf \
	-C hg38_CpG.gz
	
/sibcb2/bioinformatics/software/Miniconda3/bin/python run_snpCov.py \
	-I /sibcb1/bioinformatics/shijiantao/WGBS_Normal/alignment/Z0000043W/Z0000043W.bam \
	-S out/dbsnp_146.hg38.common_nonCpG.vcf \
	-C hg38_CpG.gz > res/Z0000043W_cov.txt
