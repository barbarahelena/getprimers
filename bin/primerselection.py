#!/usr/bin/env python3

import argparse
import pandas as pd

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Select primers based on quality and specificity criteria.")
    parser.add_argument("--primer_file", required=True, help="Path to the primer3 results file.")
    parser.add_argument("--quality_file", required=True, help="Path to the primer quality results file.")
    parser.add_argument("--specificity_file", required=True, help="Path to the primer specificity results file.")
    parser.add_argument("--output_file", required=True, help="Path to save the selected primers.")
    parser.add_argument("--meta_id", required=False, help="Metadata ID for tracking (optional).")
    args = parser.parse_args()

    # Load the primer3 results
    primer3_results = pd.read_csv(args.primer_file)

    # Load the primer quality results
    primer_quality = pd.read_csv(args.quality_file)

    # Load the off-target results
    specificity_results = pd.read_csv(args.specificity_file)

    # Merge all tables into one based on Gene and Primer_Set
    combined_results = primer3_results.merge(
        primer_quality,
        on=["Gene", "Primer_Set", "Forward_Primer", "Reverse_Primer", "Forward_Primer_Adapter", "Reverse_Primer_Adapter", "Forward_Tm", "Reverse_Tm"],
        how="inner"
    ).merge(
        specificity_results,
        on=["Gene", "Primer_Set"],
        how="inner"
    )

    # Save the combined results to a CSV file
    combined_results_file = args.output_file.replace(".csv", "_combined.csv")
    combined_results.to_csv(combined_results_file, index=False)
    print(f"Combined results saved to '{combined_results_file}'")

    # Filter primers that meet the criteria: no hairpins, no dimers, and no off-targets
    filtered_results = combined_results[
        (combined_results["Hairpin_Flag"] == "No") &
        (combined_results["Dimer_Flag"] == "No") &
        (combined_results["Has_OffTarget_Binding"] == False)
    ]

    # Save the filtered results to the specified output file
    filtered_results.to_csv(args.output_file, index=False)
    print(f"Filtered results saved to '{args.output_file}'")

if __name__ == "__main__":
    main()