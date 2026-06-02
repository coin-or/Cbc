#!/bin/bash
# Validate feasibility of all VRP solutions
# Solves instances and checks constraints, bounds, integrality, objective

echo "Validating VRP Solutions..."
echo "=================================================================="
echo ""

solutions_dir="solutions"
mkdir -p "$solutions_dir"

passed=0
failed=0
failed_instances=""

for f in cvrp_*.mps vrppd_*.mps; do
    if [ ! -f "$f" ]; then continue; fi

    name=$(basename $f .mps)
    sol_file="$solutions_dir/${name}.sol"

    printf "%-30s " "$name"

    # Solve and save solution
    result=$(timeout 120 ../src/mipster $f -sec 60 -solve -solu "$sol_file" 2>&1)

    if [ ! -f "$sol_file" ]; then
        printf "✗ NO SOLUTION\n"
        ((failed++))
        failed_instances="$failed_instances $name"
        continue
    fi

    # Validate solution feasibility
    validation=$(mipster_validate_sol $f "$sol_file" 2>&1)
    exit_code=$?

    if [ $exit_code -eq 0 ]; then
        # Extract objective from solution file
        obj=$(head -1 "$sol_file" | awk '{print $2}')
        printf "✓ FEASIBLE  obj=%s\n" "$obj"
        ((passed++))
    else
        printf "✗ INFEASIBLE\n"
        echo "$validation" | grep -E "VIOLATION|ERROR" | head -3
        ((failed++))
        failed_instances="$failed_instances $name"
    fi
done

echo ""
echo "=================================================================="
echo "Summary: $passed feasible, $failed infeasible/failed"

if [ $failed -gt 0 ]; then
    echo ""
    echo "Failed instances:$failed_instances"
    exit 1
fi
