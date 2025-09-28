#!/bin/bash

#SBATCH --job-name=clustalo
#SBATCH --time=06:00:00
#SBATCH --mem=16G
#SBATCG --cpus-per-task=4

# Print the task id.
module load Clustal-Omega/1.2.4

# Script name
SCRIPT_NAME=$(basename "$0")
echo "script name is ${SCRIPT_NAME}"
# Function to display help message
show_help() {
    cat << EOF
Usage: ${SCRIPT_NAME} [OPTIONS]



OPTIONS:
    -t, --txt_path TXTPATH    Specify the path to text file listing genes names for clustal-omega processing (required)
    -d, --dir_path DIRPATH    Specify the path to gene fasta file directory (required)
    -o, --clustalo_out ClustaloOUT    Specify the path to clustal-omega output fasta file directory 
    -h, --help                 Display this help message and exit

EXAMPLES:
    ${SCRIPT_NAME} -t clustalo/test_data/test_data_genes_with_errors.txt -d clustalo/test_data -o clustalo/test_result

EOF
}

# Initialize variables
TXTPATH=""
DIRPATH=""
ClustaloOUT=""

# Parse command line arguments using getopt
TEMP=$(getopt -o t:d:o:h --long text_path:,dir_path:,clustalo_out:,help \
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
        -t|--txt_path)
            TXTPATH="$2"
            shift 2
            ;;
        -d|--dir_path)
            DIRPATH="$2"
            shift 2
            ;;
	-o|--clustalo_out)
	    ClustaloOUT="$2"
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
if [ -z "${TXTPATH}" ]; then
    echo "Error: text path is required. Use -t or --txt_path to specify." >&2
    echo "Use $SCRIPT_NAME -h for help." >&2
    exit 1
fi

if [ -z "${DIRPATH}" ]; then
    echo "Error: directory path is required. Use -d or --dir_path to specify." >&2
    echo "Use ${SCRIPT_NAME} -h for help." >&2
    exit 1
fi

if [ -z "${ClustaloOUT}" ]; then
    echo "Error: KaKs output path is required. Use -o or --clustalo_out to specify." >&2
    echo "Use ${SCRIPT_NAME} -h for help." >&2
    exit 1
fi

# Print confirmation messages
echo "=================================="
echo ""
echo "The following arguments were entered:"
echo "  text path: ${TXTPATH}"
echo "  directory path: ${DIRPATH}"
echo "  Clustal-omega output directory path: ${ClustaloOUT}"
echo ""
echo "=================================="


FILE_LIST=${TXTPATH}
FILE_DIR=${DIRPATH}
echo "${DIRPATH}"
mkdir -p ${ClustaloOUT}
# Get the file for this array task
FILE_TO_PROCESS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${FILE_LIST}")
echo "file to process: ${FILE_TO_PROCESS}"
FILE_FULL_PATH="${FILE_DIR}/${FILE_TO_PROCESS}.fa"
RESULT_FULL_PATH="${ClustaloOUT}/${FILE_TO_PROCESS}.clustal.fa"


if [[ -z "${FILE_TO_PROCESS}" ]]; then
    echo "No file found for task ID ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID} processing: ${FILE_TO_PROCESS}"
echo "File full path ${FILE_FULL_PATH}"

base_name=$(basename ${FILE_FULL_PATH})

echo "base_name ${base_name}"
#echo JOBID: $SLURM_JOBID TASK ${SLURM_ARRAY_TASK_ID}
#echo ARRAYJOBID: ${SLURM_ARRAY_JOB_ID} ${SLURM_ARRAY_TASK_ID}
clustalo --infile ${FILE_FULL_PATH} --dealign --out ${RESULT_FULL_PATH} \
	--outfmt fa --threads 4 --seqtype dna --force






