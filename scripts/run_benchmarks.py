import subprocess
import time
import os
import glob
import json
import argparse

def get_msa_files(directory):
    files = glob.glob(os.path.join(directory, "*.msa"))
    files.sort(key=lambda x: os.path.getsize(x))
    return files

def run_command(cmd):
    start = time.time()
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    end = time.time()
    return end - start, result

def main():
    parser = argparse.ArgumentParser(description="Run TN93 benchmarks.")
    parser.add_argument("--max-time", type=float, default=10.0, help="Maximum runtime in seconds before stopping (default: 10.0)")
    args_cmd = parser.parse_args()

    msa_dir = "/Users/sergei/Development/hivtrace-secure-server/test/realistic-national-2025/all/"
    ref_bin = "/usr/local/bin/tn93"
    curr_bin = "./build/tn93"

    files = get_msa_files(msa_dir)

    print(f"{'File':<30} | {'Seqs':<8} | {'Ref (s)':<10} | {'New (s)':<10} | {'New+H (s)':<10} | {'N/R':<8} | {'NH/R':<8} | {'NH/N':<8} | {'Result':<8}")
    print("-" * 125)

    for f in files:
        fname = os.path.basename(f)
        ref_csv = f"ref_{fname}.csv"
        new_csv = f"new_{fname}.csv"
        new_h_csv = f"new_h_{fname}.csv"

        opts = "-t 0.015 -a RYSMBK -g 0.05"

        # 1. Reference run
        ref_cmd = f"{ref_bin} {opts} -o {ref_csv} {f}"
        ref_time, ref_res = run_command(ref_cmd)

        # 2. New implementation (No Hamming)
        new_cmd = f"{curr_bin} {opts} -o {new_csv} {f}"
        new_time, new_res = run_command(new_cmd)

        # 3. New implementation (With Hamming)
        new_h_cmd = f"{curr_bin} {opts} -H -o {new_h_csv} {f}"
        new_h_time, new_h_res = run_command(new_h_cmd)

        # Parse JSON for sequence count
        seq_count = "N/A"
        try:
            combined = ref_res.stdout + ref_res.stderr
            start_idx = combined.find('{')
            end_idx = combined.rfind('}')
            if start_idx != -1 and end_idx != -1:
                json_data = json.loads(combined[start_idx:end_idx+1])
                seq_count = json_data.get("Sequences", "N/A")
        except:
            pass

        # Compare result identity (Ref vs New+H)
        compare_cmd = f"python3 scripts/compare_csv.py {ref_csv} {new_h_csv}"
        comp_res = subprocess.run(compare_cmd, shell=True, capture_output=True, text=True)

        result_str = "SUCCESS" if comp_res.returncode == 0 else "FAILED"

        s1 = ref_time / new_time if new_time > 0 else 0.0
        s2 = ref_time / new_h_time if new_h_time > 0 else 0.0
        s3 = new_time / new_h_time if new_h_time > 0 else 0.0

        print(f"{fname:<30} | {str(seq_count):>8} | {ref_time:>10.3f} | {new_time:>10.3f} | {new_h_time:>10.3f} | {s1:>7.2f}x | {s2:>7.2f}x | {s3:>7.2f}x | {result_str}")
        if os.path.exists(ref_csv): os.remove(ref_csv)
        if os.path.exists(new_csv): os.remove(new_csv)
        if os.path.exists(new_h_csv): os.remove(new_h_csv)

        if ref_time > args_cmd.max_time or new_time > args_cmd.max_time or new_h_time > args_cmd.max_time:
            print(f"Stopping benchmarks as runtime exceeded {args_cmd.max_time} seconds.")
            break


if __name__ == "__main__":
    main()

