input="fastq/TK305_S6_LALL_R1_001.fastq.gz"
libid="TK305"

bowtie2_ref="/mnt/data/Sequencing/EMBL_Mouse_chimera_SNP/SNP_masked_ref/C57BL-6NJ_CAST-EiJ_Nmask.bowtie2"
cutsite_path="/mnt/data/Resources/mm10.r68/recognition_sites/mm10.r68.MboI.bed.gz"

threads=10

# FASTQ quality control
mkdir fastqc
fastqc $input -o fastqc --nogroup

# Extract flags and filter by UMI quality
mkdir fastq_hq
fbarber flag extract \
	"$input" fastq_hq/$libid.hq.fastq.gz \
	--filter-qual-output fastq_hq/$libid.lq.fastq.gz \
	--unmatched-output fastq_hq/$libid.unmatched.fastq.gz \
	--log-file fastq_hq/$libid.log \
	--pattern '^(?<umi>.{8})(?<bc>.{8})(?<cs>.{4})' \
	--flagstats bc cs --filter-qual-flags umi,30,.2 \
	--threads $threads --chunk-size 200000

# Filter by prefix
mkdir fastq_prefix
fbarber flag regex \
	fastq_hq/$libid.hq.fastq.gz fastq_prefix/$libid.fastq.gz \
	--unmatched-output fastq_prefix/$libid.unmatched.fastq.gz \
	--log-file fastq_prefix/$libid.log \
	--pattern "bc,^(?<bc>GTCGTATC){s<2}$" "cs,^(?<cs>GATC){s<2}$" \
	--threads $threads --chunk-size 200000

# Align
mkdir mapping
bowtie2 \
	-x "$bowtie2_ref" fastq_prefix/$libid.fastq.gz \
	--very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder -p $threads \
	-S mapping/$libid.sam &> mapping/$libid.log

# Filter alignment
sambamba view -S mapping/$libid.sam -f bam -t $threads > mapping/$libid.bam
rm -i mapping/$libid.sam
sambamba view mapping/$libid.bam -f bam \
	-F "mapping_quality < 30" -c -t $threads \
	> mapping/$libid.lq_count.txt
sambamba view mapping/$libid.bam -f bam -t $threads \
	-F "mapping_quality >= 30 and not secondary_alignment and not unmapped and not chimeric" \
	> mapping/$libid.clean.bam
sambamba view mapping/$libid.clean.bam -f bam -c -t $threads > mapping/$libid.clean_count.txt

# Correct aligned position
mkdir atcs
sambamba view -t $threads -h -f bam -F "reverse_strand" \
	mapping/$libid.clean.bam -o atcs/$libid.clean.revs.bam
sambamba view -t $threads atcs/$libid.clean.revs.bam | \
	convert2bed --input=sam --keep-header - > atcs/$libid.clean.revs.bed
cut -f 1-4 atcs/$libid.clean.revs.bed | tr "~" $'\t' | cut -f 1,3,7,16 | gzip \
	> atcs/$libid.clean.revs.umi.txt.gz
rm atcs/$libid.clean.revs.bam atcs/$libid.clean.revs.bed

sambamba view -t $threads -h -f bam -F "not reverse_strand" \
	mapping/$libid.clean.bam -o atcs/$libid.clean.plus.bam
sambamba view -t $threads atcs/$libid.clean.plus.bam | \
	convert2bed --input=sam --keep-header - > atcs/$libid.clean.plus.bed
cut -f 1-4 atcs/$libid.clean.plus.bed | tr "~" $'\t' | cut -f 1,2,7,16 | gzip \
	> atcs/$libid.clean.plus.umi.txt.gz
rm atcs/$libid.clean.plus.bam atcs/$libid.clean.plus.bed

# Group UMIs
scripts/group_umis.py \
	atcs/$libid.clean.plus.umi.txt.gz \
	atcs/$libid.clean.revs.umi.txt.gz \
	atcs/$libid.clean.umis.txt.gz \
	--compress-level 6 --len 4
rm atcs/$libid.clean.plus.umi.txt.gz atcs/$libid.clean.revs.umi.txt.gz

# Assign UMIs to cutsites
scripts/umis2cutsite.py \
	atcs/$libid.clean.umis.txt.gz $cutsite_path \
	atcs/$libid.clean.umis_at_cs.txt.gz --compress --threads $threads
rm atcs/$libid.clean.umis.txt.gz

# Deduplicate
mkdir dedup
scripts/umi_dedupl.R \
	atcs/$libid.clean.umis_at_cs.txt.gz \
	dedup/$libid.clean.umis_dedupd.txt.gz \
	-c 20 -r 10000

# Generate final bed
mkdir bed
zcat dedup/$libid.clean.umis_dedupd.txt.gz | \
	awk 'BEGIN{FS=OFS="\t"}{print $1 FS $2 FS $2 FS "pos_"NR FS $4}' | \
	gzip > bed > $libid.bed.gz

