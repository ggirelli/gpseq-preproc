libid="TK311"

snp_file="/mnt/data/Sequencing/EMBL_Mouse_chimera_SNP/SNPs/C57BL-6J_CAST-EiJ.txt.gz"

threads=10

SNPsplit --weird --no_sort --snp_file "$snp_file" \
	mapping/$libid.clean.bam

mkdir -p genome1/mapping genome2/mapping
mv mapping/$libid.clean.genome1.bam genome1/mapping/$libid.bam
mv mapping/$libid.clean.genome2.bam genome2/mapping/$libid.bam

# By genome -------------------------------------------------------------------------------------

glab_list=("genome1")
glab_list+=("genome2")

for glab in ${glab_list[@]}; do
	cutsite_path="/mnt/data/Sequencing/EMBL_Mouse_chimera_SNP/C57BL-6NJ_CAST-EiJ_"$glab".GG30Oct2020.MboI.sorted.bed.gz"
	
	# Correct aligned position
	mkdir -p $glab/atcs
	sambamba view -q -t $threads -h -f bam -F "reverse_strand" \
		$glab/mapping/$libid.bam -o $glab/atcs/$libid.clean.revs.bam
	sambamba view -q -t $threads $glab/atcs/$libid.clean.revs.bam | \
		convert2bed --input=sam --keep-header - > $glab/atcs/$libid.clean.revs.bed
	cut -f 1-4 $glab/atcs/$libid.clean.revs.bed | tr "~" $'\t' | cut -f 1,3,7,16 | gzip \
		> $glab/atcs/$libid.clean.revs.umi.txt.gz
	rm $glab/atcs/$libid.clean.revs.bam $glab/atcs/$libid.clean.revs.bed

	sambamba view -q -t $threads -h -f bam -F "not reverse_strand" \
		$glab/mapping/$libid.bam -o $glab/atcs/$libid.clean.plus.bam
	sambamba view -q -t $threads $glab/atcs/$libid.clean.plus.bam | \
		convert2bed --input=sam --keep-header - > $glab/atcs/$libid.clean.plus.bed
	cut -f 1-4 $glab/atcs/$libid.clean.plus.bed | tr "~" $'\t' | cut -f 1,2,7,16 | gzip \
		> $glab/atcs/$libid.clean.plus.umi.txt.gz
	rm $glab/atcs/$libid.clean.plus.bam $glab/atcs/$libid.clean.plus.bed

	# Group UMIs
	scripts/group_umis.py \
		$glab/atcs/$libid.clean.plus.umi.txt.gz \
		$glab/atcs/$libid.clean.revs.umi.txt.gz \
		$glab/atcs/$libid.clean.umis.txt.gz \
		--compress-level 6 --len 4
	rm $glab/atcs/$libid.clean.plus.umi.txt.gz $glab/atcs/$libid.clean.revs.umi.txt.gz

	# Assign UMIs to cutsites
	cutsite_path="/mnt/data/Resources/mm10.r68/recognition_sites/mm10.r68.MboI.bed.gz"
	scripts/umis2cutsite.py \
		$glab/atcs/$libid.clean.umis.txt.gz $cutsite_path \
		$glab/atcs/$libid.clean.umis_at_cs.txt.gz --compress --threads $threads
	rm $glab/atcs/$libid.clean.umis.txt.gz

	# Deduplicate
	mkdir -p $glab/dedup
	scripts/umi_dedupl.R \
		$glab/atcs/$libid.clean.umis_at_cs.txt.gz \
		$glab/dedup/$libid.clean.umis_dedupd.txt.gz \
		-c $threads -r 10000

	# Generate final bed
	mkdir -p $glab/bed
	zcat $glab/dedup/$libid.clean.umis_dedupd.txt.gz | \
		awk 'BEGIN{FS=OFS="\t"}{print $1 FS $2 FS $2 FS "pos_"NR FS $4}' | \
		gzip > $glab/bed/$libid.$glab.bed.gz
done
