#!/bin/bash
# Example benchmark runs for PandaDock using PDBbind dataset

# Set your PDBbind dataset path here
PDBBIND_DIR="/path/to/pdbbind"
OUTPUT_BASE="benchmark_results"

echo "PandaDock PDBbind Benchmark Examples"
echo "====================================="

# Check if PDBbind directory exists
if [ ! -d "$PDBBIND_DIR" ]; then
    echo "Error: PDBbind directory not found at $PDBBIND_DIR"
    echo "Please edit this script and set the correct PDBBIND_DIR path"
    exit 1
fi

# Create timestamp for unique output directories
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

echo "Starting benchmarks at $(date)"
echo "Output will be saved to ${OUTPUT_BASE}_${TIMESTAMP}_*"
echo ""

# Example 1: Quick test with 10 entries (CPU only)
echo "Example 1: Quick test with 10 entries (CPU only)"
echo "================================================"
python pdbbind_benchmark.py \
    --pdbbind_dir "$PDBBIND_DIR" \
    --output "${OUTPUT_BASE}_${TIMESTAMP}_quick_test" \
    --max_entries 10 \
    --algorithms genetic \
    --devices CPU \
    --scoring standard \
    --iterations 25 \
    --verbose

echo ""

# Example 2: CPU vs GPU performance comparison
echo "Example 2: CPU vs GPU performance comparison (50 entries)"
echo "========================================================="
python pdbbind_benchmark.py \
    --pdbbind_dir "$PDBBIND_DIR" \
    --output "${OUTPUT_BASE}_${TIMESTAMP}_cpu_vs_gpu" \
    --max_entries 50 \
    --algorithms genetic \
    --devices CPU GPU \
    --scoring standard \
    --iterations 50 \
    --parallel_jobs 2

echo ""

# Example 3: Algorithm comparison
echo "Example 3: Algorithm comparison (100 entries)"
echo "=============================================="
python pdbbind_benchmark.py \
    --pdbbind_dir "$PDBBIND_DIR" \
    --output "${OUTPUT_BASE}_${TIMESTAMP}_algorithms" \
    --max_entries 100 \
    --algorithms genetic random \
    --devices CPU \
    --scoring standard \
    --iterations 50 \
    --parallel_jobs 4

echo ""

# Example 4: Scoring function comparison
echo "Example 4: Scoring function comparison (75 entries)"
echo "==================================================="
python pdbbind_benchmark.py \
    --pdbbind_dir "$PDBBIND_DIR" \
    --output "${OUTPUT_BASE}_${TIMESTAMP}_scoring" \
    --max_entries 75 \
    --algorithms genetic \
    --devices CPU \
    --scoring standard enhanced \
    --iterations 50 \
    --parallel_jobs 3

echo ""

# Example 5: Comprehensive benchmark (smaller subset)
echo "Example 5: Comprehensive benchmark (200 entries)"
echo "================================================="
python pdbbind_benchmark.py \
    --pdbbind_dir "$PDBBIND_DIR" \
    --output "${OUTPUT_BASE}_${TIMESTAMP}_comprehensive" \
    --max_entries 200 \
    --algorithms genetic random \
    --devices CPU GPU \
    --scoring standard enhanced \
    --iterations 50 \
    --parallel_jobs 6

echo ""

# Example 6: Full scale benchmark (use with caution!)
echo "Example 6: Full scale benchmark (all entries - WARNING: This will take a long time!)"
echo "===================================================================================="
echo "Uncomment the following lines to run full benchmark:"
echo "# python pdbbind_benchmark.py \\"
echo "#     --pdbbind_dir \"$PDBBIND_DIR\" \\"
echo "#     --output \"${OUTPUT_BASE}_${TIMESTAMP}_full\" \\"
echo "#     --algorithms genetic random pandadock \\"
echo "#     --devices CPU GPU \\"
echo "#     --scoring standard enhanced physics-based \\"
echo "#     --iterations 100 \\"
echo "#     --parallel_jobs 8"

echo ""
echo "Benchmarks completed at $(date)"
echo ""
echo "Results Analysis:"
echo "=================="
echo "1. Check the output directories for detailed results"
echo "2. Look for benchmark_analysis.json for statistical summaries"
echo "3. Check the plots/ subdirectory for visualizations"
echo "4. CSV files can be imported into Excel/R/Python for further analysis"