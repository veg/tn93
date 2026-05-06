import subprocess
import json
import csv
import random
import sys
import os
import argparse

def run_tn93(binary, args, input_file, output_csv):
    cmd = [binary] + args.split()
    if "-o" not in args:
        cmd += ["-o", output_csv]
    cmd += [input_file]
    
    # Run and capture both stdout and stderr
    result = subprocess.run(cmd, capture_output=True, text=True)
    return result.stdout, result.stderr, result.returncode

def extract_subset(full_csv, subset_csv, n=1500, delimiter=','):
    try:
        with open(full_csv, 'r') as f:
            reader = csv.reader(f, delimiter=delimiter)
            header = next(reader)
            rows = list(reader)
        
        if len(rows) > n:
            subset_rows = random.sample(rows, n)
        else:
            subset_rows = rows
            
        with open(subset_csv, 'w') as f:
            writer = csv.writer(f, delimiter=delimiter)
            writer.writerow(header)
            writer.writerows(subset_rows)
    except Exception as e:
        print(f"Error extracting subset: {e}")
        sys.exit(1)

def filter_json(data, ignore_keys):
    if not ignore_keys:
        return data
    
    filtered = data.copy()
    keys_to_del = []
    for k in ignore_keys:
        for key in filtered.keys():
            if key.startswith(k):
                keys_to_del.append(key)
    
    for key in set(keys_to_del):
        filtered.pop(key, None)
    return filtered

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary", required=True)
    parser.add_argument("--ref-binary", help="Path to reference binary (required only if generating ground truth)")
    parser.add_argument("--input", required=True)
    parser.add_argument("--args", default="")
    parser.add_argument("--name", required=True)
    parser.add_argument("--ground-truth-dir", default="tests/ground_truth")
    
    args = parser.parse_args()
    
    # Extract delimiter if present in args
    delimiter = ','
    if "-D" in args.args:
        # Simplistic parsing of -D
        parts = args.args.split()
        if "-D" in parts:
            idx = parts.index("-D")
            if idx + 1 < len(parts):
                delimiter = parts[idx+1]

    if not os.path.exists(args.ground_truth_dir):
        os.makedirs(args.ground_truth_dir)
        
    ref_json_path = os.path.join(args.ground_truth_dir, f"{args.name}.json")
    ref_csv_path = os.path.join(args.ground_truth_dir, f"{args.name}.csv")
    
    if not os.path.exists(ref_json_path) or not os.path.exists(ref_csv_path):
        if not args.ref_binary:
            print(f"Error: Ground truth missing for {args.name} and no --ref-binary provided.")
            sys.exit(1)
            
        print(f"Generating ground truth for {args.name}...")
        temp_csv = f"{args.name}_full_ref.csv"
        # Strip -H from reference arguments since reference binary might not support it
        ref_args = args.args.replace("-H", "").strip()
        stdout, stderr, rc = run_tn93(args.ref_binary, ref_args, args.input, temp_csv)
        if rc != 0:
            print(f"Reference binary failed with code {rc}")
            print(stderr)
            if os.path.exists(temp_csv): os.remove(temp_csv)
            sys.exit(1)
        
        # Strip potential non-JSON lines from combined output
        combined = stdout + stderr
        json_start = combined.find('{')
        json_end = combined.rfind('}')
        if json_start != -1 and json_end != -1:
            stderr_json = combined[json_start:json_end+1]
        else:
            print("Could not find JSON in reference output")
            print("STDOUT:", stdout)
            print("STDERR:", stderr)
            sys.exit(1)

        with open(ref_json_path, 'w') as f:
            f.write(stderr_json)
        
        extract_subset(temp_csv, ref_csv_path, delimiter=delimiter)
        if os.path.exists(temp_csv): os.remove(temp_csv)

    # Run target binary
    print(f"Running test {args.name}...")
    target_csv = f"{args.name}_target.csv"
    target_stdout, target_stderr, target_rc = run_tn93(args.binary, args.args, args.input, target_csv)
    
    if target_rc != 0:
        print(f"Target binary failed with code {target_rc}")
        print(target_stderr)
        if os.path.exists(target_csv): os.remove(target_csv)
        sys.exit(1)
        
    # Compare JSON
    with open(ref_json_path, 'r') as f:
        ref_data = json.load(f)
        
    combined_target = target_stdout + target_stderr
    json_start = combined_target.find('{')
    json_end = combined_target.rfind('}')
    if json_start != -1 and json_end != -1:
        target_data = json.loads(combined_target[json_start:json_end+1])
    else:
        print("Could not find JSON in target output")
        print("STDOUT:", target_stdout)
        print("STDERR:", target_stderr)
        sys.exit(1)
        
    ignore_keys = []
    if "-H" in args.args:
        ignore_keys = ["Histogram", "Maximum distance", "Mean distance", "Note"]
    
    filtered_ref = filter_json(ref_data, ignore_keys)
    filtered_target = filter_json(target_data, ignore_keys)
    
    if filtered_ref != filtered_target:
        print("JSON output mismatch!")
        print("Reference (filtered):", json.dumps(filtered_ref, indent=2))
        print("Target (filtered):", json.dumps(filtered_target, indent=2))
        sys.exit(1)
    
    # Compare CSV subset
    with open(ref_csv_path, 'r') as f:
        reader = csv.reader(f, delimiter=delimiter)
        header = next(reader)
        # Use a dict for fast lookup. Key is (id1, id2), value is distance.
        ref_lookup = {tuple(sorted(row[:2])): row[2] for row in reader}
        
    with open(target_csv, 'r') as f:
        reader = csv.reader(f, delimiter=delimiter)
        header = next(reader)
        matches = 0
        for row in reader:
            pair = tuple(sorted(row[:2]))
            if pair in ref_lookup:
                # Compare distances as floats to handle precision differences
                try:
                    if abs(float(row[2]) - float(ref_lookup[pair])) > 1e-7:
                        print(f"Distance mismatch for {pair}: ref {ref_lookup[pair]}, target {row[2]}")
                        sys.exit(1)
                    matches += 1
                except ValueError:
                    if row[2] != ref_lookup[pair]:
                         print(f"Distance mismatch for {pair}: ref {ref_lookup[pair]}, target {row[2]}")
                         sys.exit(1)
                    matches += 1

    if matches < len(ref_lookup):
        print(f"Could not find all reference pairs in target output. Found {matches}/{len(ref_lookup)}")
        sys.exit(1)
        
    print(f"Test {args.name} passed!")
    if os.path.exists(target_csv):
        os.remove(target_csv)

if __name__ == "__main__":
    main()
