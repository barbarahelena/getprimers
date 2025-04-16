#!/usr/bin/env python3
import argparse
import os
import pandas as pd
from Bio.Blast.Applications import NcbiblastnCommandline

def parse_query_id(query_id):
    """
    Parse a primer query ID to extract gene name, primer set number, and primer type.
    
    Args:
        query_id (str): The query ID string (e.g., "GENE_Set1_F")
        
    Returns:
        tuple: (gene, primer_set, primer_type)
    """
    if not isinstance(query_id, str) or pd.isna(query_id):
        return ("", "", "")
        
    try:
        # For IDs like "GENE_Set1_F"
        if "_Set" in query_id:
            parts = query_id.split("_Set")
            gene = parts[0]
            
            # Extract set number and primer type
            remaining = parts[1]
            if "_" in remaining:
                set_type = remaining.split("_")
                primer_set = set_type[0]
                primer_type = set_type[1]
            else:
                # Handle case where format is different
                primer_set = remaining
                primer_type = ""
        else:
            # Handle alternative formats
            parts = query_id.split("_")
            if len(parts) >= 3:
                gene = parts[0]
                primer_set = parts[1]
                primer_type = parts[2]
            else:
                gene = query_id
                primer_set = ""
                primer_type = ""
                
        return (gene, primer_set, primer_type)
    except Exception:
        # Return empty values if parsing fails
        return ("", "", "")

def filter_blast_results(blast_file_or_df, primer_info=None, min_identity=80, min_length=16, max_mismatches=5):
    """
    Filter BLAST results to identify significant matches and potential off-target binding.
    
    Args:
        blast_file_or_df: Either a path to a BLAST results file (str) or a DataFrame with BLAST results
        primer_info (pd.DataFrame, optional): DataFrame containing primer information
        min_identity (int): Minimum percent identity threshold 
        min_length (int): Minimum alignment length
        max_mismatches (int): Maximum number of mismatches allowed
        
    Returns:
        pd.DataFrame: Filtered DataFrame with additional analysis columns
    """
    # Check if input is a string (file path) or DataFrame
    if isinstance(blast_file_or_df, str):  # Fixed: using parameter name properly
        # It's a file path, read it
        columns = [
            "Query_ID", "Subject_ID", "Percent_Identity", "Alignment_Length",
            "Mismatches", "Gap_Opens", "Query_Start", "Query_End",
            "Subject_Start", "Subject_End", "E_Value", "Bit_Score"
        ]
        try:
            blast_df = pd.read_csv(blast_file_or_df, sep='\t', names=columns, engine='python')
        except Exception as e:
            print(f"Error reading BLAST file: {e}")
            return pd.DataFrame()
    else:
        # It's already a DataFrame
        blast_df = blast_file_or_df
    
    # Now proceed with filtering
    filtered_results = blast_df.copy()
    
    # Filter out based on thresholds
    filtered_results = filtered_results[
        (filtered_results["Mismatches"] < max_mismatches) &
        (filtered_results["Percent_Identity"] >= min_identity) &
        (filtered_results["Alignment_Length"] >= min_length)
    ]

    # Analyze mismatches in the last 5 bases at the 3' end
    def has_3prime_mismatches(row, window=5):
        """Check for mismatches in the 3' end of the primer."""
        qseq = row.get("Query_Sequence", row.get("qseq", ""))
        sseq = row.get("Subject_Sequence", row.get("sseq", ""))
        
        if not isinstance(qseq, str) or not isinstance(sseq, str):
            return True
            
        q_3prime = qseq[-window:]
        s_3prime = sseq[-window:]
        
        mismatches = sum(1 for a, b in zip(q_3prime, s_3prime) if a != b)
        return mismatches > 0
    
    # Add 3' mismatch information if sequence data is available
    seq_cols = ['Query_Sequence', 'Subject_Sequence', 'qseq', 'sseq']
    if any(col in filtered_results.columns for col in seq_cols):
        filtered_results["Has_3prime_Mismatches"] = filtered_results.apply(has_3prime_mismatches, axis=1)
    
    return filtered_results

def find_expected_targets(blast_results, primers_df):
    """
    Find expected targets for each primer in the BLAST results.
    
    Args:
        blast_results (pd.DataFrame): DataFrame containing BLAST results
        primers_df (pd.DataFrame): DataFrame containing primer information
        
    Returns:
        pd.DataFrame: DataFrame with information about expected targets
    """
    # Initialize a DataFrame to store expected targets
    expected_targets = pd.DataFrame({
        'Gene': pd.Series(dtype='str'),
        'Primer_Set': pd.Series(dtype='str'),
        'Forward_Hit_Count': pd.Series(dtype='int'), 
        'Reverse_Hit_Count': pd.Series(dtype='int'),
        'Expected_Target': pd.Series(dtype='str'),
        'Target_Found': pd.Series(dtype='bool'),
        'Target_Identity': pd.Series(dtype='float')
    })
    
    # Group the blast results by Gene and Primer_Set
    grouped = blast_results.groupby(['Gene', 'Primer_Set'])
    
    # Process each primer set
    for (gene, primer_set), group in grouped:
        # Filter for perfect matches (100% identity)
        perfect_matches = group[group['Percent_Identity'] == 100.0]
        
        # Count forward and reverse primer hits
        forward_hits = perfect_matches[perfect_matches['Primer_Type'] == 'F'].shape[0]
        reverse_hits = perfect_matches[perfect_matches['Primer_Type'] == 'R'].shape[0]
        
        # Determine if target is found (both forward and reverse primers have perfect matches)
        target_found = forward_hits > 0 and reverse_hits > 0
        
        # Find the target with highest identity if any
        if not group.empty:
            target_identity = group['Percent_Identity'].max()
            expected_target = group.loc[group['Percent_Identity'].idxmax(), 'Subject_ID']
        else:
            target_identity = 0.0
            expected_target = "No target found"
        
        # Add to results
        expected_targets = pd.concat([expected_targets, pd.DataFrame({
            'Gene': [gene],
            'Primer_Set': [primer_set],
            'Forward_Hit_Count': [forward_hits],
            'Reverse_Hit_Count': [reverse_hits],
            'Expected_Target': [expected_target],
            'Target_Found': [target_found],
            'Target_Identity': [target_identity]
        })], ignore_index=True)
    
    return expected_targets

def find_potential_amplicons(blast_results, max_amplicon_size=1000):
    """
    Find potential amplicons from filtered BLAST results by identifying pairs of forward and reverse primers.
    
    Args:
        blast_results (pd.DataFrame): DataFrame containing filtered BLAST results
        max_amplicon_size (int): Maximum allowed amplicon size in base pairs
        
    Returns:
        pd.DataFrame: DataFrame containing potential amplicons
    """
    # Initialize an empty DataFrame to store amplicons
    amplicons = pd.DataFrame(columns=[
        'Genome', 'Subject_ID', 'Forward_Primer', 'Reverse_Primer', 
        'Amplicon_Size', 'Forward_Identity', 'Reverse_Identity', 'Is_Off_Target'
    ])
    
    # Ensure we have the required columns
    required_cols = ['Genome', 'Subject_ID', 'Query_ID', 'Gene', 'Primer_Set', 'Primer_Type',
                     'Subject_Start', 'Subject_End', 'Percent_Identity']
    if not all(col in blast_results.columns for col in required_cols):
        print("Missing required columns in BLAST results for amplicon finding")
        return amplicons
    
    # Group by Genome, Subject_ID, Gene and Primer_Set to find primer pairs
    grouped = blast_results.groupby(['Genome', 'Subject_ID', 'Gene', 'Primer_Set'])
    
    # Process each group to find forward-reverse primer pairs
    for (genome, subject_id, gene, primer_set), group in grouped:
        if len(group) < 2:
            continue  # Need at least 2 primers (forward and reverse) for an amplicon
            
        # Find forward and reverse primers
        forward_primers = group[group['Primer_Type'] == 'F']
        reverse_primers = group[group['Primer_Type'] == 'R']
        
        if forward_primers.empty or reverse_primers.empty:
            continue  # Need both forward and reverse primers
            
        # For each forward-reverse pair, calculate potential amplicon
        for _, f_primer in forward_primers.iterrows():
            for _, r_primer in reverse_primers.iterrows():
                # Calculate amplicon size
                if f_primer['Subject_Start'] < r_primer['Subject_Start']:
                    # Forward primer is upstream of reverse primer (normal orientation)
                    amplicon_size = r_primer['Subject_End'] - f_primer['Subject_Start'] + 1
                    orientation = '+'
                else:
                    # Reverse primer is upstream of forward primer (could be valid for circular DNA)
                    amplicon_size = f_primer['Subject_End'] - r_primer['Subject_Start'] + 1
                    orientation = '-'
                
                # Skip if amplicon size exceeds maximum
                if amplicon_size > max_amplicon_size or amplicon_size <= 0:
                    continue
                
                # Determine if this is an off-target amplicon more carefully:
                # An amplicon is off-target if:
                # 1. It's found on a sequence other than the target gene's expected target
                # 2. Or if either primer has less than perfect matching at the 3' end
                # 3. Or if the overall identity is below 95% (but at least 80% to be detected at all)
                
                is_off_target = False
                expected_target_key = f"{f_primer['Gene']}_{f_primer['Primer_Set']}"
                
                # Check if this subject_id is in the expected targets list for this primer set
                if expected_target_key in expected_targets_lookup:
                    expected_targets_for_set = expected_targets_lookup[expected_target_key]
                    expected_subject_ids = [t['Target'] for t in expected_targets_for_set]
                    
                    # If not in expected targets, mark as off-target
                    if subject_id not in expected_subject_ids:
                        is_off_target = True
                else:
                    # Can't determine expected target, be conservative and mark as off-target
                    is_off_target = True
                
                # If either primer has mismatches at the 3' end, mark as off-target
                if ('Has_3prime_Mismatches' in f_primer and f_primer['Has_3prime_Mismatches']) or \
                   ('Has_3prime_Mismatches' in r_primer and r_primer['Has_3prime_Mismatches']):
                    is_off_target = True
                
                # If either primer has less than 95% identity, mark as off-target
                if f_primer['Percent_Identity'] < 95 or r_primer['Percent_Identity'] < 95:
                    is_off_target = True
                
                # Add to amplicons DataFrame
                new_amplicon = pd.DataFrame({
                    'Genome': [genome],
                    'Subject_ID': [subject_id],
                    'Forward_Primer': [f_primer['Query_ID']],
                    'Reverse_Primer': [r_primer['Query_ID']],
                    'Amplicon_Size': [amplicon_size],
                    'Orientation': [orientation],
                    'Forward_Identity': [f_primer['Percent_Identity']],
                    'Reverse_Identity': [r_primer['Percent_Identity']],
                    'Is_Off_Target': [is_off_target],
                    'Gene': [gene],  # Keep track of the gene for better reporting
                    'Primer_Set': [primer_set]  # Keep track of the primer set
                })
                
                amplicons = pd.concat([amplicons, new_amplicon], ignore_index=True)
    
    # Sort by Genome and amplicon size
    if not amplicons.empty:
        amplicons = amplicons.sort_values(['Genome', 'Amplicon_Size'])
    
    return amplicons
    
    # Sort by Genome and amplicon size
    if not amplicons.empty:
        amplicons = amplicons.sort_values(['Genome', 'Amplicon_Size'])
    
    return amplicons

def summarize_amplicons(amplicons_df, max_examples=5):
    """
    Summarize amplicon data into a concise format for reporting
    
    Args:
        amplicons_df (pd.DataFrame): DataFrame containing amplicon information
        max_examples (int): Maximum number of examples to include in the summary
    
    Returns:
        str: A summary string describing the amplicons
    """
    if amplicons_df.empty:
        return "No amplicons found"
    
    # Count amplicons by size ranges
    size_ranges = {
        "small (<100bp)": (0, 100),
        "medium (100-500bp)": (100, 500),
        "large (>500bp)": (500, float('inf'))
    }
    
    summary_parts = []
    total_amplicons = len(amplicons_df)
    
    # Count by size range
    for range_name, (min_size, max_size) in size_ranges.items():
        count = amplicons_df[(amplicons_df['Amplicon_Size'] >= min_size) & 
                            (amplicons_df['Amplicon_Size'] < max_size)].shape[0]
        if count > 0:
            summary_parts.append(f"{count} {range_name}")
    
    # Add unique genome count
    unique_genomes = amplicons_df['Genome'].nunique()
    summary_parts.append(f"across {unique_genomes} genomes")
    
    # Add examples (but avoid duplicates)
    examples_df = amplicons_df.drop_duplicates(subset=['Genome', 'Subject_ID', 'Amplicon_Size'])
    if not examples_df.empty:
        examples = examples_df.sort_values('Amplicon_Size').head(max_examples)
        example_texts = []
        
        for _, row in examples.iterrows():
            example_texts.append(
                f"{row['Amplicon_Size']}bp on {row['Subject_ID']}"
            )
        
        if len(examples) < len(examples_df):
            example_texts.append(f"and {len(examples_df) - len(examples)} more")
            
        summary_parts.append("Examples: " + ", ".join(example_texts))
    
    return f"{total_amplicons} total amplicons ({'; '.join(summary_parts)})"

# Use f-strings more effectively
def get_primer_id(gene, primer_set, primer_type):
    return f"{gene}_Set{primer_set}_{primer_type}"

def main(primer_file, reference_genomes, meta_id):
    """
    Main function to perform primer specificity analysis.
    """
    # Read primer sequences from the input file
    primers = pd.read_csv(primer_file)

    # Debug: Print the columns in the primer file
    print("Columns in the primer file:", primers.columns)

    # Check if required columns exist
    required_columns = ['Gene', 'Primer_Set', 'Forward_Primer', 'Reverse_Primer']
    for col in required_columns:
        if col not in primers.columns:
            raise ValueError(f"Missing required column '{col}' in the primer file.")

    # Write primers to a FASTA file for BLAST
    primer_fasta = f"{meta_id}_primers.fasta"
    with open(primer_fasta, "w") as f:
        for index, row in primers.iterrows():
            fwd_seq = row['Forward_Primer']
            rev_seq = row['Reverse_Primer']
            
            # Write FASTA entries
            f.write(f">{row['Gene']}_Set{row['Primer_Set']}_F\n{fwd_seq}\n")
            f.write(f">{row['Gene']}_Set{row['Primer_Set']}_R\n{rev_seq}\n")

    # Create a dictionary to store primer lengths for later use
    primer_lengths = {}
    for _, row in primers.iterrows():
        f_key = get_primer_id(row['Gene'], row['Primer_Set'], 'F')
        r_key = get_primer_id(row['Gene'], row['Primer_Set'], 'R')
        primer_lengths[f_key] = len(row['Forward_Primer'])
        primer_lengths[r_key] = len(row['Reverse_Primer'])

    # Initialize a DataFrame to store combined results
    combined_results = pd.DataFrame()

    # Iterate over each reference genome and run BLAST
    reference_genomes_list = reference_genomes.split()
    for ref_genome in reference_genomes_list:
        genome_name = os.path.basename(ref_genome).replace(".fna", "")
        blast_output = f"{genome_name}_blast_results.tsv"

        # Calculate the word size based on the primer length
        primer_length = primers["Forward_Primer"].str.len().max()  # Get the maximum primer length
        word_size = min(11, max(7, int(primer_length * 0.4)))  # 40% of primer length, between 7-11

        # Construct the BLAST command
        blastn_cline = NcbiblastnCommandline(
            query=primer_fasta,
            subject=ref_genome,
            outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
            out=blast_output,
            task="blastn-short",
            word_size=word_size,
            evalue=30,  
            dust="no",           # Don't filter low-complexity regions
            soft_masking="true", # Better handling of repetitive regions
            penalty=-2,          # Standard penalty for mismatches
            reward=1,            # Standard reward for matches
            gapopen=5,           # Cost to open a gap
            gapextend=2,         # Cost to extend a gap
            max_target_seqs=5000 # Report many hits per primer
        )

        # Debug: Print the BLAST command
        print(f"Running BLAST for {genome_name} with word_size={word_size}: {blastn_cline}")

        # Execute the BLAST command
        stdout, stderr = blastn_cline()

        # Check if there are errors
        if stderr:
            print(f"BLAST stderr: {stderr}")

        # Parse BLAST results and add genome name
        columns = [
            "Query_ID", "Subject_ID", "Percent_Identity", "Alignment_Length",
            "Mismatches", "Gap_Openings", "Query_Start", "Query_End",
            "Subject_Start", "Subject_End", "E-value", "Bit_Score"
        ]
        
        try:
            # Explicitly specify engine='python' to avoid the warnings
            blast_results = pd.read_csv(blast_output, sep="\t", names=columns, engine="python")
            
            # Add genome name to the results
            blast_results["Genome"] = genome_name
            
            # Extract Gene and Primer_Set more efficiently  
            blast_results[['Gene', 'temp']] = blast_results['Query_ID'].str.split('_Set', n=1, expand=True)
            blast_results['Primer_Set'] = blast_results['temp'].str.extract(r'(\d+)')
            blast_results['Primer_Type'] = blast_results['temp'].str.extract(r'\d+_([FR])')
            blast_results.drop(columns=['temp'], inplace=True)
            
            # Append to combined results
            combined_results = pd.concat([combined_results, blast_results], ignore_index=True)
        except pd.errors.EmptyDataError:
            print(f"No BLAST hits found for {genome_name}")
        except Exception as e:
            print(f"Error reading BLAST results for {genome_name}: {e}")

    # Process the combined results to analyze specificity
    if not combined_results.empty:
        # Extract Gene, Primer_Set, and Primer_Type from Query_ID if needed
        if 'Gene' not in combined_results.columns or combined_results['Gene'].isna().any():
            extracted_info = combined_results['Query_ID'].apply(parse_query_id)
            combined_results[['Gene', 'Primer_Set', 'Primer_Type']] = pd.DataFrame(
                extracted_info.tolist(), index=combined_results.index
            )
        
        # Filter the combined results for specificity analysis
        print("Filtering BLAST results to focus on high-quality hits...")
        filtered_results = filter_blast_results(combined_results, primers)
        
        # Save filtered results for inspection (optional - can be commented out to save space)
        filtered_results.to_csv(f"{meta_id}_filtered_blast_results.tsv", sep="\t", index=False)
        print(f"  - Found {len(filtered_results)} high-quality primer binding sites")
        
        # Find expected targets (likely the gene itself) for each primer
        print("Identifying expected targets for each primer set...")
        try:
            expected_targets = find_expected_targets(filtered_results, primers)
            expected_targets.to_csv(f"{meta_id}_expected_targets.csv", index=False)
            print(f"  - Identified {expected_targets['Target_Found'].sum()} primer sets with confirmed targets")
        except Exception as e:
            print(f"Warning: Error in find_expected_targets: {e}")
            # Create an empty DataFrame as a fallback
            expected_targets = pd.DataFrame(columns=['Gene', 'Primer_Set', 'Forward_Hit_Count', 
                                                   'Reverse_Hit_Count', 'Expected_Target', 
                                                   'Target_Found', 'Target_Identity'])
        
        # Flag 3' end matches efficiently
        filtered_results['is_3prime_match'] = False
        if 'Has_3prime_Mismatches' in filtered_results.columns:
            filtered_results.loc[filtered_results['Has_3prime_Mismatches'] == False, 'is_3prime_match'] = True
        
        # Find potential amplicons using the filtered results
        print("Finding potential amplicons from high-quality primer binding sites...")
        try:
            potential_amplicons = find_potential_amplicons(filtered_results, max_amplicon_size=1000)
            print(f"  - Found {len(potential_amplicons)} potential amplicons within size limit (≤1000 bp)")
            
            # Write amplicons to file for inspection
            potential_amplicons.to_csv(f"{meta_id}_potential_amplicons.csv", index=False)
        except Exception as e:
            print(f"Warning: Error in find_potential_amplicons: {e}")
            potential_amplicons = pd.DataFrame(columns=[
                'Genome', 'Subject_ID', 'Forward_Primer', 'Reverse_Primer', 
                'Amplicon_Size', 'Forward_Identity', 'Reverse_Identity', 'Is_Off_Target'
            ])
        
        # Create a lookup dictionary for expected targets to optimize matching
        expected_targets_lookup = {}
        if not expected_targets.empty:
            for _, row in expected_targets.iterrows():
                key = f"{row['Gene']}_{row['Primer_Set']}"
                if key not in expected_targets_lookup:
                    expected_targets_lookup[key] = []
                expected_targets_lookup[key].append({
                    'Gene': row['Gene'],
                    'Primer_Set': row['Primer_Set'],
                    'Target': row['Expected_Target']
                })
        
        # Find concerning matches more efficiently
        concerning_matches = []
        high_quality_matches = filtered_results[
            (filtered_results['Percent_Identity'] >= 85) & 
            (filtered_results['is_3prime_match'] == True)
        ]
        
        if not high_quality_matches.empty:
            # Group by Query_ID to process primer by primer
            for qid, matches in high_quality_matches.groupby('Query_ID'):
                # Parse gene and primer set from query ID
                gene, primer_set, _ = parse_query_id(qid)
                target_key = f"{gene}_{primer_set}"
                
                # Get expected targets for this primer
                expected_targets_for_primer = expected_targets_lookup.get(target_key, [])
                expected_target_ids = [target['Target'] for target in expected_targets_for_primer]
                
                # Find matches not in expected targets
                for _, match in matches.iterrows():
                    if match['Subject_ID'] not in expected_target_ids:
                        concerning_matches.append(match)
        
        # Find potential amplicons with maximum size restriction (1000 bp)
        try:
            # Use filtered results to find potential amplicons, not the raw BLAST results
            potential_amplicons = find_potential_amplicons(filtered_results, max_amplicon_size=1000)
        except Exception as e:
            print(f"Warning: Error in find_potential_amplicons: {e}")
            potential_amplicons = pd.DataFrame(columns=[
                'Genome', 'Subject_ID', 'Forward_Primer', 'Reverse_Primer', 
                'Amplicon_Size', 'Forward_Identity', 'Reverse_Identity', 'Is_Off_Target'
            ])
        
        # Create a more efficient specificity summary using groupby
        # Get unique gene/primer set combinations from the primer file
        specificity_summary = primers[['Gene', 'Primer_Set']].drop_duplicates().copy()
        specificity_summary['Has_OffTarget_Binding'] = False
        specificity_summary['OffTarget_Count'] = 0
        specificity_summary['Off_Target_Amplicons'] = ""
        specificity_summary['Specificity_Flag'] = "Pass"
        specificity_summary['Notes'] = ""
        
        # Process off-target amplicons
        if not potential_amplicons.empty:
            print(f"\n🔍 Found {len(potential_amplicons)} potential amplicons within size limit (≤1000 bp)")
            off_target_amplicons = potential_amplicons[potential_amplicons['Is_Off_Target'] == True]
            
            if not off_target_amplicons.empty:
                print(f"\n⚠️ Found {len(off_target_amplicons)} potential off-target amplification sites")
                
                # Save to file
                off_target_amplicons.to_csv(f"{meta_id}_off_target_amplicons.csv", index=False)
                
                # Efficiently update specificity summary
                for gene, primer_set in specificity_summary[['Gene', 'Primer_Set']].itertuples(index=False):
                    # Get amplicons for this primer set - using more efficient string operations
                    f_pattern = get_primer_id(gene, primer_set, 'F')
                    r_pattern = get_primer_id(gene, primer_set, 'R')
                    
                    # Find amplicons associated with this primer set
                    primer_amplicons = potential_amplicons[
                        potential_amplicons['Forward_Primer'].str.startswith(f_pattern) | 
                        potential_amplicons['Reverse_Primer'].str.startswith(f_pattern) |
                        potential_amplicons['Forward_Primer'].str.startswith(r_pattern) | 
                        potential_amplicons['Reverse_Primer'].str.startswith(r_pattern)
                    ]
                    
                    # Find off-target amplicons
                    primer_off_targets = primer_amplicons[primer_amplicons['Is_Off_Target'] == True]
                    
                    # Update the summary row if any off-targets found
                    if not primer_off_targets.empty:
                        mask = (specificity_summary['Gene'] == gene) & (specificity_summary['Primer_Set'] == primer_set)
                        if mask.any():
                            specificity_summary.loc[mask, 'Has_OffTarget_Binding'] = True
                            specificity_summary.loc[mask, 'OffTarget_Count'] = len(primer_off_targets)
                            
                            # Use summarize_amplicons function to create concise summary
                            try:
                                specificity_summary.loc[mask, 'Off_Target_Amplicons'] = summarize_amplicons(primer_off_targets)
                                specificity_summary.loc[mask, 'Notes'] = summarize_amplicons(primer_amplicons, max_examples=3)
                            except Exception as e:
                                print(f"Warning: Error in summarize_amplicons: {e}")
                                # Fallback to simple count
                                specificity_summary.loc[mask, 'Off_Target_Amplicons'] = f"{len(primer_off_targets)} amplicons"
                                
                            specificity_summary.loc[mask, 'Specificity_Flag'] = "Fail"
            else:
                print("\n✅ No off-target amplification detected.")
                # Create empty off_target_amplicons file
                pd.DataFrame(columns=[
                    'Genome', 'Subject_ID', 'Forward_Primer', 'Reverse_Primer', 
                    'Amplicon_Size', 'Forward_Identity', 'Reverse_Identity', 'Is_Off_Target'
                ]).to_csv(f"{meta_id}_off_target_amplicons.csv", index=False)
        else:
            print("\n✅ No potential amplicons found within size limit (≤1000 bp)")
            # Create empty off_target_amplicons file
            pd.DataFrame(columns=[
                'Genome', 'Subject_ID', 'Forward_Primer', 'Reverse_Primer', 
                'Amplicon_Size', 'Forward_Identity', 'Reverse_Identity', 'Is_Off_Target'
            ]).to_csv(f"{meta_id}_off_target_amplicons.csv", index=False)
        
        # Convert Primer_Set to string to ensure consistent types
        specificity_summary['Primer_Set'] = specificity_summary['Primer_Set'].astype(str)
        
        # Save results
        specificity_summary.to_csv(f"{meta_id}_specificity_summary.csv", index=False)
        combined_results.to_csv(f"{meta_id}_specificity_results.csv", index=False)
        
        print(f"\nSpecificity summary saved to {meta_id}_specificity_summary.csv")
        print(f"Full results saved to {meta_id}_specificity_results.csv")
        
        # Print summary statistics
        pass_count = (specificity_summary['Specificity_Flag'] == "Pass").sum()
        fail_count = (specificity_summary['Specificity_Flag'] == "Fail").sum()
        print(f"\nSummary: {pass_count} primer sets passed, {fail_count} primer sets failed specificity check")
        
        # Print amplicon statistics
        if not potential_amplicons.empty:
            on_target_count = (potential_amplicons['Is_Off_Target'] == False).sum()
            off_target_count = len(off_target_amplicons) if 'off_target_amplicons' in locals() else 0
            print(f"Amplicons: {on_target_count} on-target, {off_target_count} off-target")
    else:
        print("No BLAST hits found for any primer. Check your reference genomes or BLAST parameters.")
        # Create empty output files
        pd.DataFrame(columns=[
            'Gene', 'Primer_Set', 'Has_OffTarget_Binding', 'OffTarget_Count', 
            'Off_Target_Amplicons', 'Specificity_Flag', 'Notes'
        ]).to_csv(f"{meta_id}_specificity_summary.csv", index=False)
        
        pd.DataFrame(columns=columns + [
            'Genome', 'Gene', 'Primer_Set', 'Primer_Type', 'is_3prime_match'
        ]).to_csv(f"{meta_id}_specificity_results.csv", index=False)
        
        pd.DataFrame(columns=[
            'Genome', 'Subject_ID', 'Forward_Primer', 'Reverse_Primer', 
            'Amplicon_Size', 'Forward_Identity', 'Reverse_Identity', 'Is_Off_Target'
        ]).to_csv(f"{meta_id}_off_target_amplicons.csv", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Perform primer specificity analysis.")
    parser.add_argument("--primer_file", required=True, help="Path to the primer file.")
    parser.add_argument("--reference_genomes", required=True, help="Paths to the reference genomes (space-separated).")
    parser.add_argument("--meta_id", required=True, help="Meta ID for the analysis.")

    args = parser.parse_args()

    main(args.primer_file, args.reference_genomes, args.meta_id)