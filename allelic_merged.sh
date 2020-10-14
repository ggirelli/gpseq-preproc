libid_list=("TK300")
libid_list+=("TK301")
libid_list+=("TK302")
libid_list+=("TK303")
libid_list+=("TK304")
libid_list+=("TK305")

snp_file="/mnt/data/Sequencing/EMBL_Mouse_chimera_SNP/SNPs/C57BL-6J_CAST-EiJ.txt.gz"

threads=10

mkdir -p genome12/mapping
for libid in ${libid_list[@]}; do
	sambamba sort genome1/mapping/$libid.bam -t $threads -p -q
	sambamba sort genome2/mapping/$libid.bam -t $threads -p -q
	sambamba merge genome12/mapping/$libid.bam genome1/mapping/$libid.sorted.bam genome2/mapping/$libid.sorted.bam \
		-t $threads -p -q
	rm genome1/mapping/$libid.sorted.bam genome2/mapping/$libid.sorted.bam
	rm genome1/mapping/$libid.sorted.bam.bai genome2/mapping/$libid.sorted.bam.bai
done

for libid in ${libid_list[@]}; do
	# Correct aligned position
	mkdir genome12/atcs
	sambamba view -q -t $threads -h -f bam -F "reverse_strand" \
		genome12/mapping/$libid.bam -o genome12/atcs/$libid.clean.revs.bam
	sambamba view -q -t $threads genome12/atcs/$libid.clean.revs.bam | \
		convert2bed --input=sam --keep-header - > genome12/atcs/$libid.clean.revs.bed
	cut -f 1-4 genome12/atcs/$libid.clean.revs.bed | tr "~" $'\t' | cut -f 1,3,7,16 | gzip \
		> genome12/atcs/$libid.clean.revs.umi.txt.gz
	rm genome12/atcs/$libid.clean.revs.bam genome12/atcs/$libid.clean.revs.bed

	sambamba view -q -t $threads -h -f bam -F "not reverse_strand" \
		genome12/mapping/$libid.bam -o genome12/atcs/$libid.clean.plus.bam
	sambamba view -q -t $threads genome12/atcs/$libid.clean.plus.bam | \
		convert2bed --input=sam --keep-header - > genome12/atcs/$libid.clean.plus.bed
	cut -f 1-4 genome12/atcs/$libid.clean.plus.bed | tr "~" $'\t' | cut -f 1,2,7,16 | gzip \
		> genome12/atcs/$libid.clean.plus.umi.txt.gz
	rm genome12/atcs/$libid.clean.plus.bam genome12/atcs/$libid.clean.plus.bed

	# Group UMIs
	scripts/group_umis.py \
		genome12/atcs/$libid.clean.plus.umi.txt.gz \
		genome12/atcs/$libid.clean.revs.umi.txt.gz \
		genome12/atcs/$libid.clean.umis.txt.gz \
		--compress-level 6 --len 4
	rm genome12/atcs/$libid.clean.plus.umi.txt.gz genome12/atcs/$libid.clean.revs.umi.txt.gz

	# Assign UMIs to cutsites
	cutsite_path="/mnt/data/Resources/mm10.r68/recognition_sites/mm10.r68.MboI.bed.gz"
	scripts/umis2cutsite.py \
		genome12/atcs/$libid.clean.umis.txt.gz $cutsite_path \
		genome12/atcs/$libid.clean.umis_at_cs.txt.gz --compress --threads $threads
	rm genome12/atcs/$libid.clean.umis.txt.gz

	# Deduplicate
	mkdir genome12/dedup
	scripts/umi_dedupl.R \
		genome12/atcs/$libid.clean.umis_at_cs.txt.gz \
		genome12/dedup/$libid.clean.umis_dedupd.txt.gz \
		-c 20 -r 10000

	# Generate final bed
	mkdir genome12/bed
	zcat genome12/dedup/$libid.clean.umis_dedupd.txt.gz | \
		awk 'BEGIN{FS=OFS="\t"}{print $1 FS $2 FS $2 FS "pos_"NR FS $4}' | \
		gzip > genome12/bed/$libid.bed.gz
done
