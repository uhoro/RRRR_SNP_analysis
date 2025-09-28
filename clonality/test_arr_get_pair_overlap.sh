#!/bin/bash

#SBATCH --job-name=array
#SBATCH --time=03:00:00
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCG --cpus-per-task=4


# Script name
SCRIPT_NAME=$(basename "$0")
echo "script name is ${SCRIPT_NAME}"
# Function to display help message
show_help() {
    cat << EOF
Usage: ${SCRIPT_NAME} [OPTIONS]

This script computes clonality for given pairs of samples

OPTIONS:
    -x, --txt_path1 TXTPATH1    Specify the path to text file listing sample 1 (required)
    -y, --txt_path2 TXTPATH2    Specify the path to text file listing sample 2 (required)
    -d, --data_path DATAPATH    Specify the path to the variant data (required)
    -r, --outdir_path OUTPATH    Specify the path to the output (required)
    -h, --help                 Display this help message and exit

EXAMPLES:
    ${SCRIPT_NAME} -t -d

EOF
}

# Initialize variables
TXTPATH1=""
TXTPATH2=""
DATAPATH=""
OUTPATH=""

# Parse command line arguments using getopt
TEMP=$(getopt -o x:y:d:r:h --long txt_path1:,txt_path2:,data_path:,outdir_path:,help \
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
        -x|--txt_path1)
            TXTPATH1="$2"
            shift 2
            ;;
        -y|--txt_path2)
            TXTPATH2="$2"
            shift 2
            ;;
        -d|--data_path)
            DATAPATH="$2"
            shift 2
            ;;
        -r|--outdir_path)
            OUTPATH="$2"
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
if [ -z "${TXTPATH1}" ]; then
    echo "Error: text path for sample1 is required. Use -x or --txt_path1 to specify." >&2
    echo "Use $SCRIPT_NAME -h for help." >&2
    exit 1
fi

if [ -z "${TXTPATH2}" ]; then
    echo "Error: text path for sample2 is required. Use -y or --txt_path2 to specify." >&2
    echo "Use ${SCRIPT_NAME} -h for help." >&2
    exit 1
fi

if [ -z "${DATAPATH}" ]; then
    echo "Error: data path is required. Use -d or --data_path to specify." >&2
    echo "Use ${SCRIPT_NAME} -h for help." >&2
    exit 1
fi

if [ -z "${OUTPATH}" ]; then
    echo "Error: output directory path is required. Use -r or --outdir_path to specify." >&2
    echo "Use ${SCRIPT_NAME} -h for help." >&2
    exit 1
fi

WD=### USER DEFINED WD ### 

txt1full="${WD}${TXTPATH1}"
txt2full="${WD}${TXTPATH2}"
datafull="${WD}${DATAPATH}"
outdirfull="${WD}${OUTPATH}"

# Print confirmation messages
echo "=================================="
echo ""
echo "The following arguments were entered:"
echo "  sample 1 text path: ${txt1full}"
echo "  sample 2 text path: ${txt2full}"
echo "  data path: ${datafull}"
echo "  outdir path: ${outdirfull}"
echo ""
echo "=================================="

# Get the file for this array task
SAMPLE1_TO_PROCESS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${txt1full}")

# Get the file for this array task
SAMPLE2_TO_PROCESS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${txt2full}")

if [[ -z "${SAMPLE1_TO_PROCESS}" ]]; then
    echo "No sample 1 found for task ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

if [[ -z "${SAMPLE2_TO_PROCESS}" ]]; then
    echo "No sample 2 found for task ID $SLURM_ARRAY_TASK_ID"
    exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID} processing sample1 : ${SAMPLE1_TO_PROCESS}"
echo "Task ${SLURM_ARRAY_TASK_ID} processing sample2 : ${SAMPLE2_TO_PROCESS}"



eval "$(micromamba shell hook --shell bash)"
env_path=### USER DEFINED PATH TO MICROMAMBA ENVIRONMENT FOR RRRR_ANALYSIS ###
micromamba activate ${env_path}
python get_pair_overlap_250827.py --data_fname ${datafull} --s1 ${SAMPLE1_TO_PROCESS} --s2 ${SAMPLE2_TO_PROCESS} --pair_overlap_dir ${outdirfull}
micromamba deactivate


