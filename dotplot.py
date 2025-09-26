### Authors: Ragib Ashan,  and Laura Katz (helped with AI)
### Updated on january 25th of 2025 by Marie Leleu

### This script is about dotplot visualization, intially writen for the Ciliates project
# it includes both dotplot and graph of divergences in alignment to compare 2 contigs between each other
# 3 flags available: -dotplot, -window, -both

### Input: the script takes a -i flag for input folder with the contigs to dotplot in fasta format (ending with .fasta)
# and in which  there is another subfolder called for_sliding_window for the -window flag

### Output: the script creates an Output folder with all pdf and png created
# files are named like the Input files (removing the .fasta part)

### running: python compare.py -i InputFolder -both

#!/usr/bin/env python3

import os
import argparse
from Bio import SeqIO, AlignIO
import matplotlib.pyplot as plt
import numpy as np
from Bio.Seq import Seq
import glob

def setup_output_directory():
    """Create output directory if it doesn't exist"""
    if not os.path.exists('Output'):
        os.makedirs('Output')

def read_sequences(fasta_file):
    """Read sequences from FASTA file"""
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    if len(sequences) != 2:
        raise ValueError(f"Input FASTA file {fasta_file} must contain exactly two sequences")
    return sequences[0], sequences[1]

def create_dotplot(seq1, seq2, window_size=10):
    """Generate dotplot data for two sequences"""
    matches_forward = []
    matches_reverse = []
    seq2_rev = seq2.reverse_complement()
    
    for i in range(len(seq1.seq) - window_size + 1):
        window1 = str(seq1.seq[i:i+window_size])
        for j in range(len(seq2.seq) - window_size + 1):
            window2 = str(seq2.seq[j:j+window_size])
            window2_rev = str(seq2_rev.seq[j:j+window_size])
            
            if window1 == window2:
                matches_forward.append((j, i))
            if window1 == window2_rev:
                matches_reverse.append((len(seq2.seq) - j - window_size, i))
    
    return matches_forward, matches_reverse

def plot_dotplot(seq1, seq2, matches_forward, matches_reverse, output_base):
    """Create and save dotplot"""
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Plot forward matches in blue
    if matches_forward:
        x_forward, y_forward = zip(*matches_forward)
        ax.scatter(x_forward, y_forward, c='blue', s=1, alpha=0.5, label='Forward')
    
    # Plot reverse matches in red
    if matches_reverse:
        x_reverse, y_reverse = zip(*matches_reverse)
        ax.scatter(x_reverse, y_reverse, c='red', s=1, alpha=0.5, label='Reverse')
    
    # Set axis labels and ticks
    ax.set_xlabel(f'{seq2.id} (bp)', fontname='Arial')
    ax.set_ylabel(f'{seq1.id} (bp)', fontname='Arial')
    
    # Set ticks every 1000bp
    xticks = np.arange(0, len(seq2.seq), 1000)
    yticks = np.arange(0, len(seq1.seq), 1000)
    ax.set_xticks(xticks)
    ax.set_yticks(yticks)
    ax.tick_params(labelsize=8)
    
    plt.legend()
    
    # Save plots
    plt.savefig(f'Output/{output_base}_dotplot.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'Output/{output_base}_dotplot.png', format='png', dpi=300, bbox_inches='tight')
    plt.close()

def calculate_pairwise_divergence(alignment, window_size=50, step=1):
    """Calculate average pairwise divergence in sliding windows"""
    seq_count = len(alignment)
    align_length = alignment.get_alignment_length()
    divergence = []
    
    for i in range(0, align_length - window_size + 1, step):
        window_divergence = 0
        comparisons = 0
        
        # Compare each pair of sequences
        for j in range(seq_count):
            for k in range(j + 1, seq_count):
                window_j = str(alignment[j].seq[i:i+window_size])
                window_k = str(alignment[k].seq[i:i+window_size])
                
                # Count differences, excluding gaps
                diff_count = sum(1 for x, y in zip(window_j, window_k) 
                               if x != y and x != '-' and y != '-')
                valid_pos = sum(1 for x, y in zip(window_j, window_k) 
                              if x != '-' and y != '-')
                
                if valid_pos > 0:
                    window_divergence += diff_count / valid_pos
                    comparisons += 1
        
        if comparisons > 0:
            divergence.append(window_divergence / comparisons)
        else:
            divergence.append(0)
    
    return divergence

def plot_window_divergence(alignment_file, window_size=50, step=1):
    """Create and save window divergence plot for an alignment"""
    # Read alignment
    alignment = AlignIO.read(alignment_file, "fasta")
    
    # Calculate divergence
    divergence = calculate_pairwise_divergence(alignment, window_size, step)
    
    # Create plot
    fig, ax = plt.subplots(figsize=(12, 6))
    
    x = range(0, len(divergence) * step, step)
    ax.plot(x, divergence, 'b-')
    
    ax.set_xlabel('Position (bp)', fontname='Arial')
    ax.set_ylabel('Average pairwise divergence', fontname='Arial')
    ax.set_title(os.path.basename(alignment_file), fontname='Arial')
    
    # Save plots
    base_name = os.path.splitext(os.path.basename(alignment_file))[0]
    plt.savefig(f'Output/{base_name}_divergence.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig(f'Output/{base_name}_divergence.png', format='png', dpi=300, bbox_inches='tight')
    plt.close()

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Compare sequences using dotplot and window similarity')
    parser.add_argument('-i', '--input_dir', required=True, help='Input directory containing FASTA files')
    parser.add_argument('-dotplot', action='store_true', help='Generate dotplot only')
    parser.add_argument('-window', action='store_true', help='Generate window similarity plot only')
    parser.add_argument('-both', action='store_true', help='Generate both plots')
    
    args = parser.parse_args()
    
    # Create output directory
    setup_output_directory()
    
    # Get all FASTA files from input directory
    fasta_files = glob.glob(os.path.join(args.input_dir, '*.fasta'))
    if not fasta_files:
        print(f"No FASTA files found in directory: {args.input_dir}")
        return
    
    # Perform requested analyses
    if args.dotplot or args.both:
        for fasta_file in fasta_files:
            try:
                base_name = os.path.splitext(os.path.basename(fasta_file))[0]
                seq1, seq2 = read_sequences(fasta_file)
                matches_forward, matches_reverse = create_dotplot(seq1, seq2)
                plot_dotplot(seq1, seq2, matches_forward, matches_reverse, base_name)
                print(f"Dotplot analysis complete for {fasta_file}")
            except ValueError as e:
                print(f"Error processing {fasta_file}: {str(e)}")
                continue
    
    if args.window or args.both:
        window_dir = os.path.join(args.input_dir, 'for_sliding_window')
        if os.path.exists(window_dir):
            alignment_files = glob.glob(os.path.join(window_dir, '*.fasta'))
            if not alignment_files:
                print("No alignment files found in for_sliding_window directory")
            else:
                for alignment_file in alignment_files:
                    try:
                        plot_window_divergence(alignment_file)
                        print(f"Window divergence analysis complete for {alignment_file}")
                    except Exception as e:
                        print(f"Error processing {alignment_file}: {str(e)}")
                        continue
        else:
            print(f"Directory not found: {window_dir}")
    
    if not (args.dotplot or args.window or args.both):
        print("Please specify analysis type: -dotplot, -window, or -both")

if __name__ == "__main__":
    main()
