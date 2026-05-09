import subprocess
import time
import os
import glob
import json

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
    msa_dir = "/Users/sergei/Development/hivtrace-secure-server/test/realistic-national-2025/all/"
    ref_bin = "/usr/local/bin/tn93"
    curr_bin = "./build/tn93"
    
    files = get_msa_files(msa_dir)
    
    print(f"{'File':<30} | {'Seqs':<8} | {'Ref (s)':<10} | {'Curr (s)':<10} | {'Speedup':<8} | {'Result':<8}")
    print("-" * 90)
    
    for f in files:
        fname = os.path.basename(f)
        ref_csv = f"ref_{fname}.csv"
        curr_csv = f"curr_{fname}.csv"
        
        opts = "-t 0.015 -a RYSMBK -g 0.05"
        ref_cmd = f"{ref_bin} {opts} -o {ref_csv} {f}"
        ref_time, ref_res = run_command(ref_cmd)
        
        curr_cmd = f"{curr_bin} {opts} -H -o {curr_csv} {f}"
        curr_time, curr_res = run_command(curr_cmd)

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
        
        compare_cmd = f"python3 scripts/compare_csv.py {ref_csv} {curr_csv}"
        comp_res = subprocess.run(compare_cmd, shell=True, capture_output=True, text=True)
        
        result_str = "SUCCESS" if comp_res.returncode == 0 else "FAILED"
        speedup = ref_time / curr_time if curr_time > 0 else 0.0
        
        print(f"{fname:<30} | {str(seq_count):>8} | {ref_time:>10.3f} | {curr_time:>10.3f} | {speedup:>7.2f}x | {result_str}")
        
        if os.path.exists(ref_csv): os.remove(ref_csv)
        if os.path.exists(curr_csv): os.remove(curr_csv)
        
        if ref_time > 120 or curr_time > 120:
            print("Stopping benchmarks as runtime exceeded 2 minutes.")
            break

if __name__ == "__main__":
    main()
