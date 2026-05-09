import csv
import sys
import argparse

def load_csv(filename, delimiter=','):
    data = {}
    try:
        with open(filename, 'r') as f:
            reader = csv.reader(f, delimiter=delimiter)
            header = next(reader)
            for row in reader:
                if len(row) < 3:
                    continue
                pair = tuple(sorted([row[0], row[1]]))
                distance = float(row[2])
                
                if pair in data:
                    data[pair].append(distance)
                else:
                    data[pair] = [distance]
                    
        for pair in data:
            data[pair].sort()
            
        return data, header
    except Exception as e:
        print(f"Error loading {filename}: {e}")
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description="Compare two TN93 CSV files for identity up to row permutation.")
    parser.add_argument("file1", help="First CSV file")
    parser.add_argument("file2", help="Second CSV file")
    parser.add_argument("-d", "--delimiter", default=",", help="CSV delimiter (default: ,)")
    parser.add_argument("-t", "--tolerance", type=float, default=1e-9, help="Floating point tolerance for distances (default: 1e-9)")
    
    args = parser.parse_args()
    
    data1, header1 = load_csv(args.file1, args.delimiter)
    data2, header2 = load_csv(args.file2, args.delimiter)
    
    if len(data1) != len(data2):
        print(f"FAILED: Row count mismatch. {args.file1} has {len(data1)} unique pairs, {args.file2} has {len(data2)}.")
        sys.exit(1)
        
    for pair, distances1 in data1.items():
        if pair not in data2:
            print(f"FAILED: Pair {pair} found in {args.file1} but not in {args.file2}.")
            sys.exit(1)
            
        distances2 = data2[pair]
        if len(distances1) != len(distances2):
             print(f"FAILED: Pair {pair} has different number of occurrences. {len(distances1)} vs {len(distances2)}.")
             sys.exit(1)
             
        for d1, d2 in zip(distances1, distances2):
            if abs(d1 - d2) > args.tolerance:
                print(f"FAILED: Distance mismatch for pair {pair}: {d1} vs {d2} (diff: {abs(d1-d2)})")
                sys.exit(1)
                
    print("SUCCESS: Files are identical up to row permutation.")

if __name__ == "__main__":
    main()
