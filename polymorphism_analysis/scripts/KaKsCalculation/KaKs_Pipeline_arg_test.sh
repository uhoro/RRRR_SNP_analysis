#!/bin/bash
#SBATCH --job-name=pipeline_master

# Script name
SCRIPT_NAME=$(basename "$0")
echo "script name is ${SCRIPT_NAME}"
# Function to display help message
show_help() {
    cat << EOF
Usage: ${SCRIPT_NAME} [OPTIONS]

This script calculates KaKs values from variant axt files 

OPTIONS:
    -a, --axt_dir_path AXTDIR    Specify the path to axt file directory (required)
    -i, --int_dir_path INTDIR    Specify the path to the intermediate output directory (required)
    -o, --kaks_dir KaKsDIR    Specify the path to kaks file directory (required)
    -x, --sample1 SAMPLE1	Specify the sample 1 id (required)
    -y, --sample2 SAMPLE2	Specify the sample 2 id (required)
    -h, --help                 Display this help message and exit

EXAMPLES:

sbatch -J pipe -c 4 -e KaKs_Pipeline_arg_test.err -o KaKs_Pipeline_arg_test.out KaKs_Pipeline_arg_test.sh -a KaKs_pipeline_data_250622/axt_files_strict_250622 -i KaKs_pipeline_out_250622 -o KaKs_pipeline_result_250622/strict -x C1G11 -y C1C5v2
    

EOF
}

# Initialize variables
AXTDIR=""
INTDIR=""
KaKsDIR=""
SAMPLE1=""
SAMPLE2=""

# Parse command line arguments using getopt
TEMP=$(getopt -o a:i:o:x:y:h --long axt_dir_path:,int_dir_path:,kaks_dir:,sample1:,sample2:,help \
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
        -a|--axt_dir_path)
            AXTDIR="$2"
            shift 2
            ;;
        -i|--int_dir_path)
            INTDIR="$2"
            shift 2
            ;;
	-o|--kaks_dir)
	    KaKsDIR="$2"
	    shift 2
	    ;;
	-x|--sample1)
	    SAMPLE1="$2"
	    shift 2
	    ;;
	-y|--sample2)
	    SAMPLE2="$2"
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
if [ -z "${AXTDIR}" ]; then
    echo "Error: axt dir is required. Use -a or --axt_dir_path to specify." >&2
    echo "Use $SCRIPT_NAME -h for help." >&2
    exit 1
fi

if [ -z "${INTDIR}" ]; then
    echo "Error: directory path is required. Use -i or --int_dir_path to specify." >&2
    echo "Use ${SCRIPT_NAME} -h for help." >&2
    exit 1
fi

if [ -z "${KaKsDIR}" ]; then
    echo "Error: KaKs output path is required. Use -o or --kaks_dir to specify." >&2
    echo "Use ${SCRIPT_NAME} -h for help." >&2
    exit 1
fi

if [ -z "${SAMPLE1}" ]; then
    echo "Error: sample1 id is required. Use -x or --sample1 to specify." >&2
    echo "Use ${SCRIPT_NAME} -h for help." >&2
    exit 1
fi

if [ -z "${SAMPLE2}" ]; then
    echo "Error: sample2 id is required. Use -y or --sample2 to specify." >&2
    echo "Use ${SCRIPT_NAME} -h for help." >&2
    exit 1
fi

# Print confirmation messages
echo "=================================="
echo ""
echo "The following arguments were entered:"
echo "  axt dir path: ${AXTDIR}"
echo "  intermediate directory path: ${INTDIR}"
echo "  kaks result directory path: ${KaKsDIR}"
echo "  sample1 id : ${SAMPLE1}"
echo "  sample2 id : ${SAMPLE2}"
echo ""
echo "=================================="


# Variables

SCRIPT_DIR= ### User defined directory ###
AXT_DIR_PATH="${SCRIPT_DIR}/${AXTDIR}"
#AXT_DIR_PATH="${SCRIPT_DIR}/KaKs_pipeline_data_250615/rg_axt_files"
WD_inter_OUT="${SCRIPT_DIR}/${INTDIR}"
#WD_inter_OUT="${SCRIPT_DIR}/KaKs_pipeline_out_250615"
mkdir -p ${WD_inter_OUT}
KaKs_DIR="${SCRIPT_DIR}/${KaKsDIR}"
#KaKs_DIR="${SCRIPT_DIR}/KaKs_pipeline_result_250615_rg"
#SAMPLE1="C1G11"
#SAMPLE2="C1C5v2"
TotalFileNumber=14000

# Specify directory for java kakscalculator2 function
WDSPLIT=### User defined directory ###



SAMPLE_PAIR="${SAMPLE1}_${SAMPLE2}"
sample_axt_path="${AXT_DIR_PATH}/${SAMPLE_PAIR}"
sample_axt_list_path="${AXT_DIR_PATH}/${SAMPLE_PAIR}_axt_list.txt"
sample_KaKs_OUT="${KaKs_DIR}/${SAMPLE_PAIR}"
mkdir -p ${sample_KaKs_OUT}


shopt -s extglob
echo -n "" > ${sample_axt_list_path}
for f in ${sample_axt_path}/*; do
        if [[ $f =~ ^.*_FUN_[0-9]{6}\.axt$ ]]; then
                base_name=$(basename $f)
                echo ${base_name} >> ${sample_axt_list_path}
        fi
done


linenumber=$(wc -l < ${sample_axt_list_path})
if [[ ${linenumber} -ge 1 ]]; then 
        echo "yes! Suitable axt file found in ${sample_axt_path}" 
else 
     	echo "no suitable axt file found in ${sample_axt_path}. Exit"
	exit 6 
fi 


# Get split axt file
JOB1=$(sbatch --parsable --chdir=${WDSPLIT} --wait --array=1-${TotalFileNumber}%20 -c 20 \
		--output="${WD_inter_OUT}/array_%A_%a.out" \
		-e="${WD_inter_OUT}/array_%A_%a.err" \
		 test_SplitArray_arg.sh -t ${sample_axt_list_path} -d ${sample_axt_path})

#wait # it's not waiting for JOB1 to finish...

# Get split axt file list
sample_split_axt_list_path="${AXT_DIR_PATH}/${SAMPLE_PAIR}_split_axt_list.txt"
#ls ${sample_axt_path} | grep -E '*_FUN_[0-9]{6}split_114_6\.axt' > ${sample_split_axt_list_path}
shopt -s extglob
echo -n "" > ${sample_split_axt_list_path}
for f in ${sample_axt_path}/*; do
        if [[ $f =~ ^.*_FUN_[0-9]{6}split_114_6\.axt$ ]]; then
                base_name=$(basename $f)
                echo ${base_name} >> ${sample_split_axt_list_path}
        fi
done

linenumber=$(wc -l < ${sample_split_axt_list_path})
if [[ ${linenumber} -ge 1 ]]; then 
        echo "split axt files found" 
else 
     	echo "no suitable split axt file found in ${sample_axt_path}. Exit"
        exit 7
fi 



echo "sample split axt list created: ${sample_split_axt_list_path}"
echo "start the following job sbatch --parsable --dependency=afterany:${JOB1} --array=1-30%10 \
		--output="${WD_inter_OUT}/kaks_%A_%a.out" \
		-e="${WD_inter_OUT}/kaks_%A_%a.err" \
		"${SCRIPT_DIR}/test_KaKsArray_arg.sh" -t ${sample_split_axt_list_path} -d ${sample_axt_path} \
		-o ${sample_KaKs_OUT}"

# Get kaks analysis results for split axt files
JOB2=$(sbatch --parsable --dependency=afterany:${JOB1} --wait --array=1-${TotalFileNumber}%20 -c 20 \
		--output="${WD_inter_OUT}/kaks_%A_%a.out" \
		-e="${WD_inter_OUT}/kaks_%A_%a.err" \
		"${SCRIPT_DIR}/test_KaKsArray_arg.sh" -t ${sample_split_axt_list_path} -d ${sample_axt_path} \
		-o ${sample_KaKs_OUT})








