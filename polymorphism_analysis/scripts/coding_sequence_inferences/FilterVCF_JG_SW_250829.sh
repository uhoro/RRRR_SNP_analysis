#!/bin/bash

module load java/1.8.0
module load GenomeAnalysisTK/4.1.9.0
module load samtools/1.18
module load tabix/0.2.6

#sample="C1C5v2"
#OUTDIR="FilteredVCFs_JG_SW_250829/test"


# Script name
SCRIPT_NAME=$(basename "$0")
echo "script name is ${SCRIPT_NAME}"
# Function to display help message
show_help() {
    cat << EOF
Usage: ${SCRIPT_NAME} [OPTIONS]

This script retain sample id info only and filter out missing variants and output directory of filtered vcf.

OPTIONS:
    -i, --input input    Specify input vcf file
    -s, --sample sample    Specify the sample id (required)
    -d, --dir_path OUTDIR    Specify the path to filtered vcf file directory (required)
    -h, --help                 Display this help message and exit

EXAMPLES:
    ${SCRIPT_NAME} -i FilteredVCF_JG_250827/...vcf.gz  -s C1G11 -d FilteredVCF_JG_SW_250829/test

EOF
}

# Initialize variables
sample=""
OUTDIR=""

# Parse command line arguments using getopt
TEMP=$(getopt -o i:s:d:h --long input:,sample:,dir_path:,help \
              -n "$SCRIPT_NAME" -- "$@")

# Check for getopt errors
if [ $? != 0 ]; then
    echo "Error: Invalid arguments. Use -h for help." >&2
    exit 1
fi

# Note the quotes around '$TEMP': they are essential!
eval set -- "$TEMP"

# Extract options and their arguments
while true; do
    case "$1" in
        -i|--input)
            input="$2"
            shift 2
            ;;
        -s|--sample)
            sample="$2"
            shift 2
            ;;
        -d|--dir_path)
            OUTDIR="$2"
            shift 2
            ;;
        -h|--help)
            show_help
            exit 0
            ;;
        --)
            shift
            break
            ;;
        *)
            echo "Error: Unknown option $1" >&2
            exit 1
            ;;
    esac
done

# Check if required arguments are provided
if [ -z "${input}" ]; then
    echo "Error: input is required. Use -i or --input to specify." >&2
    echo "Use $SCRIPT_NAME -h for help." >&2
    exit 1
fi

if [ -z "${sample}" ]; then
    echo "Error: sample id is required. Use -s or --sample to specify." >&2
    echo "Use $SCRIPT_NAME -h for help." >&2
    exit 1
fi

if [ -z "${OUTDIR}" ]; then
    echo "Error: directory path is required. Use -d or --dir_path to specify." >&2
    echo "Use ${SCRIPT_NAME} -h for help." >&2
    exit 1
fi

# Print confirmation messages
echo "=================================="
echo ""
echo "The following arguments were entered:"
echo "  input: ${input}"
echo "  sample: ${sample}"
echo "  directory path: ${OUTDIR}"
echo ""
echo "=================================="

mkdir -p ${OUTDIR}

#exclude_s_list="${OUTDIR}/${sample}_excl.args"
#result_s_list="${OUTDIR}/${sample}_res.txt"

#bcftools query -l ${input} | sed "/^${sample}$/d" >> ${exclude_s_list}

init_var_sites=$(gatk CountVariants -V ${input})

gatk SelectVariants -O ${OUTDIR}/${sample}_only.vcf.gz \
	-V ${input} \
	-sn ${sample}

sample_var_sites=$(gatk CountVariants -V ${OUTDIR}/${sample}_only.vcf.gz)

gatk SelectVariants -O ${OUTDIR}/${sample}_filt.vcf.gz \
	-V ${OUTDIR}/${sample}_only.vcf.gz \
	--max-nocall-number 0

bcftools sort -o ${OUTDIR}/${sample}_filt.sorted.vcf.gz \
	${OUTDIR}/${sample}_filt.vcf.gz

gatk IndexFeatureFile -I ${OUTDIR}/${sample}_filt.sorted.vcf.gz

filt_var_sites=$(gatk CountVariants -V ${OUTDIR}/${sample}_filt.vcf.gz)

bcftools query -l "${OUTDIR}/${sample}_only.vcf.gz"

bcftools view --threads 4 -e 'GT[*] = "mis"' ${OUTDIR}/${sample}_only.vcf.gz > ${OUTDIR}/${sample}_filt_bcf.vcf
bcftools view ${OUTDIR}/${sample}_filt_bcf.vcf -Oz -o ${OUTDIR}/${sample}_filt_bcf.vcf.gz
tabix ${OUTDIR}/${sample}_filt_bcf.vcf.gz

bcf_filt_var_sites=$(gatk CountVariants -V ${OUTDIR}/${sample}_filt_bcf.vcf.gz)

echo "initial variant sites: ${init_var_sites}"
echo "sample variant sites: ${sample_var_sites}"
echo "final variant sites: ${filt_var_sites}"
echo "bcf final variant sites: ${bcf_filt_var_sites}"

