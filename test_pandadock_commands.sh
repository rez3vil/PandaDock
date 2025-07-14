#!/bin/bash

# PandaDock Command Testing Script
# Tests all combinations of commands to verify outputs

set -e  # Exit on any error

# Color codes for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Test variables
PROTEIN="tests/beta-2_alpha-1.pdb"
LIGAND="tests/propofol.pdb"
CENTER="-15.7 -17.7 8.18"
SIZE="40 40 40"

# Create test results directory
TEST_DIR="test_results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$TEST_DIR"

echo -e "${YELLOW}🧪 Starting PandaDock Command Testing${NC}"
echo -e "${YELLOW}📁 Test directory: $TEST_DIR${NC}"
echo ""

# Function to run test and check results
run_test() {
    local test_name="$1"
    local cmd="$2"
    local expected_files="$3"
    
    echo -e "${YELLOW}Testing: $test_name${NC}"
    echo "Command: $cmd"
    
    # Create test subdirectory
    local test_subdir="$TEST_DIR/$test_name"
    mkdir -p "$test_subdir"
    
    # Run the command
    if eval "$cmd --out $test_subdir"; then
        echo -e "${GREEN}✅ Command executed successfully${NC}"
        
        # Check for expected files
        local all_files_exist=true
        for file in $expected_files; do
            if [[ -f "$test_subdir/$file" ]]; then
                echo -e "${GREEN}✅ Found: $file${NC}"
            else
                echo -e "${RED}❌ Missing: $file${NC}"
                all_files_exist=false
            fi
        done
        
        # Check for poses directory
        if [[ -d "$test_subdir/poses" ]]; then
            local pose_count=$(ls "$test_subdir/poses"/*.pdb 2>/dev/null | wc -l)
            echo -e "${GREEN}✅ Poses directory found with $pose_count PDB files${NC}"
            
            # Check poses summary
            if [[ -f "$test_subdir/poses/poses_summary.csv" ]]; then
                echo -e "${GREEN}✅ Found poses_summary.csv${NC}"
                # Check for identical IC50 values (the bug we're tracking)
                local unique_ic50s=$(tail -n +2 "$test_subdir/poses/poses_summary.csv" | cut -d',' -f7 | sort -u | wc -l)
                if [[ $unique_ic50s -eq 1 ]]; then
                    echo -e "${RED}⚠️  WARNING: All poses have identical IC50 values${NC}"
                else
                    echo -e "${GREEN}✅ IC50 values vary across poses${NC}"
                fi
            fi
        else
            echo -e "${RED}❌ Poses directory not found${NC}"
            all_files_exist=false
        fi
        
        # Check for save-complex functionality
        if [[ "$cmd" == *"--save-complex"* ]]; then
            local complex_files=$(find "$test_subdir" -name "*complex*" | wc -l)
            if [[ $complex_files -gt 0 ]]; then
                echo -e "${GREEN}✅ Complex files found${NC}"
            else
                echo -e "${RED}❌ WARNING: --save-complex specified but no complex files found${NC}"
            fi
        fi
        
        if $all_files_exist; then
            echo -e "${GREEN}✅ Test PASSED${NC}"
        else
            echo -e "${RED}❌ Test FAILED - Missing expected files${NC}"
        fi
    else
        echo -e "${RED}❌ Command FAILED${NC}"
    fi
    
    echo "----------------------------------------"
    echo ""
}

# Test 1: Basic docking with default settings
run_test "basic_docking" \
    "python -m pandadock --protein $PROTEIN --ligand $LIGAND --center -15.7 -17.7 8.18 --size 40 40 40" \
    "pandadock_report.html"

# Test 2: Balanced mode with ML rescoring
run_test "balanced_ml_rescoring" \
    "python -m pandadock --protein $PROTEIN --ligand $LIGAND --mode balanced --center -15.7 -17.7 8.18 --size 40 40 40 --ml-rescoring" \
    "pandadock_report.html"

# Test 3: Precise mode (physics-based)
run_test "precise_mode" \
    "python -m pandadock --protein $PROTEIN --ligand $LIGAND --mode precise --center -15.7 -17.7 8.18 --size 40 40 40" \
    "pandadock_report.html"

# Test 4: Fast mode (GA-based)
run_test "fast_mode" \
    "python -m pandadock --protein $PROTEIN --ligand $LIGAND --mode fast --center -15.7 -17.7 8.18 --size 40 40 40" \
    "pandadock_report.html"

# Test 5: JSON output format
run_test "json_output" \
    "python -m pandadock --protein $PROTEIN --ligand $LIGAND --center -15.7 -17.7 8.18 --size 40 40 40 --report-format json" \
    "pandadock_report.json"

# Test 6: CSV output format
run_test "csv_output" \
    "python -m pandadock --protein $PROTEIN --ligand $LIGAND --center -15.7 -17.7 8.18 --size 40 40 40 --report-format csv" \
    "pandadock_report.csv"

# Test 7: Save complex functionality (the issue you reported)
run_test "save_complex" \
    "python -m pandadock --protein $PROTEIN --ligand $LIGAND --center -15.7 -17.7 8.18 --size 40 40 40 --save-complex" \
    "pandadock_report.html"

# Test 8: Flexible residues
run_test "flexible_residues" \
    "python -m pandadock --protein $PROTEIN --ligand $LIGAND --center -15.7 -17.7 8.18 --size 40 40 40 --flexible-residues 'ASN265'" \
    "pandadock_report.html"

# Test 9: Different scoring functions
for scoring in "pandacore" "pandaml" "pandaphysics"; do
    run_test "scoring_$scoring" \
        "python -m pandadock --protein $PROTEIN --ligand $LIGAND --center -15.7 -17.7 8.18 --size 40 40 40 --scoring $scoring" \
        "pandadock_report.html"
done

# Test 10: Multiple poses
run_test "multiple_poses" \
    "python -m pandadock --protein $PROTEIN --ligand $LIGAND --center -15.7 -17.7 8.18 --size 40 40 40 --num-poses 20" \
    "pandadock_report.html"

# Test 11: High exhaustiveness
run_test "high_exhaustiveness" \
    "python -m pandadock --protein $PROTEIN --ligand $LIGAND --center -15.7 -17.7 8.18 --size 40 40 40 --exhaustiveness 16" \
    "pandadock_report.html"

# Test 12: Combined options (most likely to reveal bugs)
run_test "combined_options" \
    "python -m pandadock --protein $PROTEIN --ligand $LIGAND --mode balanced --center -15.7 -17.7 8.18 --size 40 40 40 --flexible-residues 'ASN265' --ml-rescoring --save-complex --num-poses 15 --report-format csv" \
    "pandadock_report.csv"

echo -e "${YELLOW}🔍 Analyzing Results...${NC}"

# Summary analysis
total_tests=0
passed_tests=0

for test_dir in "$TEST_DIR"/*/; do
    if [[ -d "$test_dir" ]]; then
        total_tests=$((total_tests + 1))
        test_name=$(basename "$test_dir")
        
        # Check if main output file exists
        if [[ -f "$test_dir/pandadock_report.html" ]] || [[ -f "$test_dir/pandadock_report.json" ]] || [[ -f "$test_dir/pandadock_report.csv" ]]; then
            passed_tests=$((passed_tests + 1))
        fi
    fi
done

echo ""
echo -e "${YELLOW}📊 Test Summary${NC}"
echo "Total tests: $total_tests"
echo "Passed tests: $passed_tests"
echo "Failed tests: $((total_tests - passed_tests))"

# Identify common issues
echo ""
echo -e "${YELLOW}🔍 Common Issues Analysis${NC}"

# Check for IC50 issues across all tests
echo "Checking for IC50 calculation issues..."
ic50_issue_count=0
for test_dir in "$TEST_DIR"/*/; do
    if [[ -f "$test_dir/poses/poses_summary.csv" ]]; then
        unique_ic50s=$(tail -n +2 "$test_dir/poses/poses_summary.csv" | cut -d',' -f7 | sort -u | wc -l)
        if [[ $unique_ic50s -eq 1 ]]; then
            ic50_issue_count=$((ic50_issue_count + 1))
        fi
    fi
done

if [[ $ic50_issue_count -gt 0 ]]; then
    echo -e "${RED}⚠️  IC50 calculation issue found in $ic50_issue_count tests${NC}"
else
    echo -e "${GREEN}✅ IC50 calculations appear to be working correctly${NC}"
fi

# Check for save-complex issues
echo "Checking for --save-complex issues..."
save_complex_issue_count=0
for test_dir in "$TEST_DIR"/*/; do
    test_name=$(basename "$test_dir")
    if [[ "$test_name" == *"save_complex"* ]] || [[ "$test_name" == *"combined_options"* ]]; then
        complex_files=$(find "$test_dir" -name "*complex*" | wc -l)
        if [[ $complex_files -eq 0 ]]; then
            save_complex_issue_count=$((save_complex_issue_count + 1))
        fi
    fi
done

if [[ $save_complex_issue_count -gt 0 ]]; then
    echo -e "${RED}⚠️  --save-complex functionality issue found in $save_complex_issue_count tests${NC}"
else
    echo -e "${GREEN}✅ --save-complex functionality appears to be working${NC}"
fi

echo ""
echo -e "${YELLOW}📁 Test results saved in: $TEST_DIR${NC}"
echo -e "${YELLOW}🔍 Review individual test directories for detailed analysis${NC}"

# Create a summary report
cat > "$TEST_DIR/test_summary.md" << EOF
# PandaDock Command Testing Summary

## Test Results
- Total tests: $total_tests
- Passed tests: $passed_tests
- Failed tests: $((total_tests - passed_tests))

## Issues Identified
- IC50 calculation issues: $ic50_issue_count tests affected
- --save-complex functionality issues: $save_complex_issue_count tests affected

## Recommendations
1. If IC50 values are identical across poses, check the ML rescoring implementation
2. If --save-complex is not working, check the complex file saving logic
3. Review the scoring function implementations for proper affinity prediction
4. Ensure all output formats (HTML, JSON, CSV) are generating correctly

## Test Directory Structure
EOF

for test_dir in "$TEST_DIR"/*/; do
    if [[ -d "$test_dir" ]]; then
        echo "- $(basename "$test_dir")" >> "$TEST_DIR/test_summary.md"
    fi
done

echo -e "${GREEN}✅ Testing completed! Check $TEST_DIR/test_summary.md for detailed results.${NC}"