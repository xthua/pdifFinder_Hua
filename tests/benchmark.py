#!/usr/bin/env python3
"""
Benchmark script for comparing performance before and after optimizations.
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

from tests.test_performance import benchmark_varying_lengths, load_seeds, create_test_sequence, run_scanning
import PdifFinder.pdifFinder as pf

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Benchmark pdifFinder scanning performance')
    parser.add_argument('--compare-before-after', action='store_true',
                        help='Compare original vs optimized performance')
    parser.add_argument('--lengths', type=int, nargs='+', default=[1000, 5000, 10000, 20000],
                        help='Sequence lengths to test')
    parser.add_argument('--output', type=str, help='Output file for results (CSV)')
    args = parser.parse_args()
    
    if args.compare_before_after:
        print("Comparing original vs optimized scanning algorithm...")
        seedList = load_seeds()
        results = []
        for L in args.lengths:
            seq = create_test_sequence(L)
            t1 = run_scanning(pf.findMatchFragmentThreadOriginal, seedList, seq)
            t2 = run_scanning(pf.findMatchFragmentThread, seedList, seq)
            speedup = t1 / t2
            results.append((L, t1, t2, speedup))
            print(f"Length {L}: original {t1:.3f}s, optimized {t2:.3f}s, speedup {speedup:.2f}x")
        
        if args.output:
            import csv
            with open(args.output, 'w', newline='') as f:
                writer = csv.writer(f)
                writer.writerow(['length', 'time_original', 'time_optimized', 'speedup'])
                writer.writerows(results)
            print(f"Results saved to {args.output}")
        
        # Calculate average speedup
        avg_speedup = sum(r[3] for r in results) / len(results)
        print(f"Average speedup: {avg_speedup:.2f}x")
        
        # Check if optimization meets 20% improvement threshold
        if avg_speedup >= 1.2:
            print("✅ Optimization meets 20% improvement target")
        else:
            print("⚠️  Optimization does not meet 20% improvement target")
            print("   Consider further optimizations")
    
    else:
        # Default behavior: run benchmark with varying lengths
        benchmark_varying_lengths()

if __name__ == '__main__':
    main()