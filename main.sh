# input="fastq/TK270_S1_LALL_R1_001.fastq.gz"
# libid="TK270"

# bowtie2_ref="/mnt/data/Sequencing/EMBL_Mouse_chimera_SNP/SNP_masked_ref/C57BL-6NJ_CAST-EiJ_Nmask.GG20Oct2020"
# cutsite_path="/mnt/data/Sequencing/EMBL_Mouse_chimera_SNP/C57BL-6NJ_CAST-EiJ_genome12.GG30Oct2020.MboI.sorted.bed.gz"

bowtie2_ref="/mnt/data/Resources/references/GRCm38.r95.dna/Mus_musculus.GRCm38.95.dna.fa"
cutsite_path="/mnt/data/Resources/mm10.r95/recognition_sites/mm10.r95.MboI.bed.gz"

threads=10

# Filter by prefix
mkdir -p fastq_prefix
fbarber flag regex \
	fastq_hq/$libid.hq.fastq.gz fastq_prefix/$libid.fastq.gz \
	--unmatched-output fastq_prefix/$libid.unmatched.fastq.gz \
	--log-file fastq_prefix/$libid.log \
	--pattern "bc,^(?<bc>GTCGTCGA){s<2}$" "cs,^(?<cs>GATC){s<2}$" \
	--threads $threads --chunk-size 200000

# Align
mkdir -p mapping
bowtie2 \
	-x "$bowtie2_ref" fastq_prefix/$libid.fastq.gz \
	--very-sensitive -L 20 --score-min L,-0.6,-0.2 --end-to-end --reorder -p $threads \
	-S mapping/$libid.sam &> mapping/$libid.mapping.log

# Filter alignment
sambamba view -q -S mapping/$libid.sam -f bam -t $threads > mapping/$libid.bam


curl -X POST -H 'Content-type: application/json' \
	--data '{"text":"Finished mapping '"$libid"'! Waiting for user input."}' \
	https://hooks.slack.com/services/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

rm -i mapping/$libid.sam
sambamba view -q mapping/$libid.bam -f bam \
	-F "mapping_quality<30" -c -t $threads \
	> mapping/$libid.lq_count.txt
sambamba view -q mapping/$libid.bam -f bam \
	-F "ref_name=='chrM'" -c -t $threads \
	> mapping/$libid.chrM.txt
sambamba view -q mapping/$libid.bam -f bam -t $threads \
	-F "mapping_quality>=30 and not secondary_alignment and not unmapped and not chimeric and ref_name!='chrM'" \
	> mapping/$libid.clean.bam
sambamba view -q mapping/$libid.clean.bam -f bam -c -t $threads > mapping/$libid.clean_count.txt

# Correct aligned position
mkdir -p atcs
sambamba view -q -t $threads -h -f bam -F "reverse_strand" \
	mapping/$libid.clean.bam -o atcs/$libid.clean.revs.bam
sambamba view -q -t $threads atcs/$libid.clean.revs.bam | \
	convert2bed --input=sam --keep-header - > atcs/$libid.clean.revs.bed
cut -f 1-4 atcs/$libid.clean.revs.bed | tr "~" $'\t' | cut -f 1,3,7,16 | gzip \
	> atcs/$libid.clean.revs.umi.txt.gz
rm atcs/$libid.clean.revs.bam atcs/$libid.clean.revs.bed

sambamba view -q -t $threads -h -f bam -F "not reverse_strand" \
	mapping/$libid.clean.bam -o atcs/$libid.clean.plus.bam
sambamba view -q -t $threads atcs/$libid.clean.plus.bam | \
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
mkdir -p dedup
scripts/umi_dedupl.R \
	atcs/$libid.clean.umis_at_cs.txt.gz \
	dedup/$libid.clean.umis_dedupd.txt.gz \
	-c $threads -r 10000

# Generate final bed
mkdir -p bed
zcat dedup/$libid.clean.umis_dedupd.txt.gz | \
	awk 'BEGIN{FS=OFS="\t"}{print $1 FS $2 FS $2 FS "pos_"NR FS $4}' | \
	gzip > bed/$libid.bed.gz


curl -X POST -H 'Content-type: application/json' \
	--data '{"text":"Finished processing '"$libid"'!"}' \
	https://hooks.slack.com/services/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
