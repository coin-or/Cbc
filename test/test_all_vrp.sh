#!/bin/bash
# Test all VRP instances with adaptive timeouts
# Goal: find at least one solution for validation

echo "Testing all VRP instances..."
echo "=================================================================="
echo ""

passed=0
failed=0

# Adaptive timeout based on instance size
get_timeout() {
    local instance=$1
    case $instance in
        *_large*|cvrp_medium|cvrp_geo_clustered|cvrp_capacity_medium|cvrp_high_utilization|cvrp_many_small)
            echo 60  # Large/hard instances: 60s
            ;;
        *_capacity_tight*|*_demand_outliers*|cvrp_asymmetric|cvrp_few_large)
            echo 30  # Medium difficulty: 30s
            ;;
        *)
            echo 15  # Easy instances: 15s
            ;;
    esac
}

for f in cvrp_*.mps vrppd_*.mps; do
    if [ ! -f "$f" ]; then continue; fi

    name=$(basename $f .mps)
    timeout_limit=$(get_timeout $name)

    printf "%-30s [%2ds] " "$name" "$timeout_limit"

    result=$(timeout $((timeout_limit + 5)) ../src/mipster $f -sec $timeout_limit -solve 2>&1)

    if echo "$result" | grep -q "Optimal solution found"; then
        obj=$(echo "$result" | grep "Objective value:" | awk '{print $3}')
        time=$(echo "$result" | grep "Total time" | tail -1 | awk '{print $5}')
        nodes=$(echo "$result" | grep "Enumerated nodes:" | awk '{print $3}')
        printf "✓ OPTIMAL  obj=%-10s time=%6ss nodes=%s\n" "$obj" "$time" "$nodes"
        ((passed++))
    elif echo "$result" | grep -q "Stopped on time limit"; then
        # Check if we found at least one solution
        if echo "$result" | grep -q "BestSol:"; then
            obj=$(echo "$result" | grep "BestSol:" | awk '{print $2}')
            bound=$(echo "$result" | grep "Bound:" | awk '{print $2}')
            gap=$(echo "$result" | grep "Gap:" | awk '{print $2}')
            printf "⏱ TIMEOUT  obj=%-10s bound=%-10s gap=%s\n" "$obj" "$bound" "$gap"
            ((passed++))  # Has solution, validation possible
        else
            printf "⏱ TIMEOUT  (no solution found)\n"
            ((failed++))
        fi
    elif echo "$result" | grep -q "INFEASIBLE"; then
        printf "✗ INFEASIBLE\n"
        ((failed++))
    else
        printf "? UNKNOWN\n"
        ((failed++))
    fi
done

echo ""
echo "=================================================================="
echo "Summary: $passed passed (solvable), $failed failed"
echo ""
echo "Note: TIMEOUT with solution counts as passed (can validate obj)"
