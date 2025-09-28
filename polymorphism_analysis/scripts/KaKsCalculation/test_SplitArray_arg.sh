#!/bin/bash

#SBATCH --job-name=array
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1G
#SBATCG --cpus-per-task=1

# Print the task id.
# module load KaKs_Calculator/3.0
module load java

# Script name
SCRIPT_NAME=$(basename "$0")
echo "script name is ${SCRIPT_NAME}"
# Function to display help message
show_help() {
    cat << EOF
Usage: ${SCRIPT_NAME} [OPTIONS]

This script takes username and password as arguments and displays them.

OPTIONS:
    -t, --txt_path TXTPATH    Specify the path to text file listing axt file names (required)
    -d, --dir_path DIRPATH    Specify the path to axt file directory (required)
    -h, --help                 Display this help message and exit

EXAMPLES:
    ${SCRIPT_NAME} -t /pasteur/zeus/projets/p02/ecb_bioinfo/uhoro/RRRR_KaKs/test_data_pipeline -d /pasteur/zeus/projets/p02/ecb_bioinfo/uhoro/RRRR_KaKs/test_data_pipeline

EOF
}

# Initialize variables
TXTPATH=""
DIRPATH=""

# Parse command line arguments using getopt
TEMP=$(getopt -o t:d:h --long text_path:,dir_path:,help \
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

# Print confirmation messages
echo "=================================="
echo ""
echo "The following arguments were entered:"
echo "  text path: ${TXTPATH}"
echo "  directory path: ${DIRPATH}"
echo ""
echo "=================================="



FILE_LIST=${TXTPATH}
FILE_DIR=${DIRPATH}

# Get the file for this array task
FILE_TO_PROCESS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILE_LIST")

if [[ -z "${FILE_TO_PROCESS}" ]]; then
    echo "No file found for task ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID} processing: ${FILE_TO_PROCESS}"
FILE_FULL_PATH="${FILE_DIR}/${FILE_TO_PROCESS}"
echo "File full path ${FILE_FULL_PATH}"

java split ${FILE_FULL_PATH} 114 6
