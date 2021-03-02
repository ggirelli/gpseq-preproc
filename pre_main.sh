input="fastq/KG35_S6_LALL_R1_001.fastq.gz"
libid="KG35"

threads=10

# FASTQ quality control
mkdir -p fastqc
fastqc $input -o fastqc --nogroup

# Extract flags and filter by UMI quality
mkdir -p fastq_hq
fbarber flag extract \
	"$input" fastq_hq/$libid.hq.fastq.gz \
	--filter-qual-output fastq_hq/$libid.lq.fastq.gz \
	--unmatched-output fastq_hq/$libid.unmatched.fastq.gz \
	--log-file fastq_hq/$libid.log \
	--pattern 'umi8bc8cs4' --simple-pattern \
	--flagstats bc cs --filter-qual-flags umi,30,.2 \
	--threads $threads --chunk-size 200000

curl -X POST -H 'Content-type: application/json' \
	--data '{"text":"Finished running pre-main for '"$libid"'"}' \
	https://hooks.slack.com/services/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
