func() {
    echo "Usage:"
    echo "add_prefix [-f bam_dir] [-i ID] [-o out_dir]"
}

while getopts ':h:f:i:o:' OPT; do
    case $OPT in
        f) bam_dir=${OPTARG};;
        i) id="$OPTARG";;
        o) out_dir=${OPTARG};;
        h) func;;
        ?) func;;
    esac
done
java -jar ~/software/picard.jar AddOrReplaceReadGroups I=${bam_dir} O=${out_dir} RGID=$id RGLB=$id RGPL=illumina RGSM=$id RGPU=unit1