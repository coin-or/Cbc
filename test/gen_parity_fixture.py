#!/usr/bin/env python3
import os

def generate_parity_lp():
    lp_content = []
    lp_content.append("Minimize")
    # Minimize sum of binary variables
    lp_content.append("  " + " + ".join(f"x{i}" for i in range(1, 13)))
    
    lp_content.append("Subject To")
    # 10 equations to meet the minimum threshold of 10 parity rows.
    # We use a lower-triangular system:
    # row1: x1 - 2 y1 = 1
    # row_k: x_{k-1} + x_k - 2 y_k = b_k
    lp_content.append("  row1: x1 - 2 y1 = 1")
    for k in range(2, 11):
        rhs = 1 if k % 2 == 1 else 0
        lp_content.append(f"  row{k}: x{k-1} + x{k} - 2 y{k} = {rhs}")
        
    # Add non-parity rows that will NOT be removed by parity presolve
    lp_content.append("  non_parity1: 3 x11 + x12 >= 1")
    lp_content.append("  non_parity2: x10 + x11 - x12 <= 1")
        
    lp_content.append("Bounds")
    for i in range(1, 13):
        lp_content.append(f"  0 <= x{i} <= 1")
    for i in range(1, 11):
        lp_content.append(f"  0 <= y{i} <= 5")
        
    lp_content.append("Binary")
    lp_content.append("  " + " ".join(f"x{i}" for i in range(1, 13)))
    
    lp_content.append("General")
    lp_content.append("  " + " ".join(f"y{i}" for i in range(1, 11)))
    
    lp_content.append("End")
    
    # Write to file
    fixture_dir = "test/fixtures"
    os.makedirs(fixture_dir, exist_ok=True)
    filepath = os.path.join(fixture_dir, "parity_test.lp")
    with open(filepath, "w") as f:
        f.write("\n".join(lp_content) + "\n")
    print(f"Generated {filepath}")

if __name__ == "__main__":
    generate_parity_lp()
