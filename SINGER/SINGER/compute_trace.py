import tskit
import numpy as np
import pandas as pd
import argparse

def mse(x, y):
    x = np.array(x)
    y = np.array(y)
    return np.mean((x - y)**2)

def count_incompatibility(ts):
    unmapped_sites = 0
    for tree in ts.trees():
        for site in tree.sites():
            num_mutations = len(site.mutations)
            if num_mutations > 2:
                unmapped_sites += 1
                print(site.position)
            elif num_mutations == 2:
                if site.mutations[0].node != tree.root and site.mutations[1].node != tree.root:
                    unmapped_sites += 1
                    print(site.position)
    return unmapped_sites

def incompatibility_trace(prefix, indices):
    counts = []
    for index in indices:
        file_name = f"{prefix}_{index}.trees"
        try:
            ts = tskit.load(file_name)
            count = count_incompatibility(ts)
            counts.append(count)
        except FileNotFoundError:
            print(f"File not found: {file_name}")
            counts.append(None)
    return counts

def diversity_fit_mse(ts, m):
    windows = np.arange(0, ts.sequence_length, 1e6)
    windows = np.append(windows, ts.sequence_length)
    site_diversity = ts.diversity(windows=windows, mode='site')
    branch_diversity = ts.diversity(windows=windows, mode='branch')*m
    fit_mse = mse(site_diversity, branch_diversity)
    return fit_mse

def diversity_fit_trace(prefix, m, indices):
    fit_mses = []
    for index in indices:
        file_name = f"{prefix}_{index}.trees"
        try:
            ts = tskit.load(file_name)
            fit_mse = diversity_fit_mse(ts, m)
            fit_mses.append(fit_mse)
        except FileNotFoundError:
            print(f"File not found: {file_name}")
            fit_mses.append(None)
    return fit_mses


def main():
    parser = argparse.ArgumentParser(description='Compute traces for MCMC samples from SINGER.')

    parser.add_argument('-prefix', type=str, required=True, help='Prefix of the ARG sample.')
    parser.add_argument('-m', type=float, help='Mutation rate.')
    parser.add_argument('-start_index', type=int, help='The start index of the ARG sample')
    parser.add_argument('-end_index', type=int, required=True, help='The end index of the ARG sample')
    parser.add_argument('-output_filename', type=str, required=True, help='Output filename of the MCMC traces.')   
 
    args = parser.parse_args()

    indices = range(args.start_index, args.end_index + 1)
    
    # Compute incompatibility trace
    incompatibility_counts = incompatibility_trace(args.prefix, indices)
    
    # Compute diversity fit trace (if mutation rate provided)
    if args.m is not None:
        diversity_fits = diversity_fit_trace(args.prefix, args.m, indices)
    
    # Build and save results
    df = pd.DataFrame({'index': list(indices), 'incompatibility': incompatibility_counts})
    if args.m is not None:
        df['diversity_fit_mse'] = diversity_fits
    
    df.to_csv(args.output_filename, index=False)
    print(f"Traces saved to {args.output_filename}")
     
     

if __name__ == "__main__":
    main()
