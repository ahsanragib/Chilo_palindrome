#Author, date: Laura, Ragib, And AI. January 21 2025. Updated on July 3, 2025
#Motivation: Dotplot and Sliding Window (to vidualize the divergence)
#Intent: Create Dotplot of the given fasta file with sequences, and create Divergence Images on the input fasta file pairwise alignments
#Dependencies: Python3, matplotlib, argparse, numpy, glob
#Inputs: a fasta file for dotplot (Input FASTA file must contain exactly two sequences)
        #Pairwise aligned fasta files inside a folder call- 'for_sliding_window'

#Outputs: A foldes (Output). Inside it saves the dotplot and sliding windows images

#Example: 
#For Dotplot Analysis:
        # Using a single FASTA file with 2 sequences
        #python3 <script> -dotplot -f path/to/your/file.fasta
#For Window Analysis:
        # Using a directory containing alignment files
        #python3 <script> -window -d path/to/alignment/directory


#!/usr/bin/env python3

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
        raise ValueError("Input FASTA file must contain exactly two sequences")
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

def plot_dotplot(seq1, seq2, matches_forward, matches_reverse):
    """Create and save dotplot"""
    fig, ax = plt.subplots(figsize=(10, 10))

    # Plot forward matches in black
    if matches_forward:
        x_forward, y_forward = zip(*matches_forward)
        ax.scatter(x_forward, y_forward, c='black', s=1, alpha=0.5, label='Forward')

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
    plt.savefig('Output/dotplot.pdf', format='pdf', dpi=300, bbox_inches='tight')
    plt.savefig('Output/dotplot.png', format='png', dpi=300, bbox_inches='tight')
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

def plot_window_divergence(alignment_files, window_size=50, step=1):
    """Create and save window divergence plots with bases separated from X-axis"""
    divergences = []
    filenames = []

    # Read and calculate divergence for each alignment
    for alignment_file in alignment_files:
        alignment = AlignIO.read(alignment_file, "fasta")
        divergence = calculate_pairwise_divergence(alignment, window_size, step)
        divergences.append(divergence)
        filenames.append(os.path.splitext(os.path.basename(alignment_file))[0])

    # Determine global Y-axis limits
    global_min = min(min(div) for div in divergences)
    global_max = max(max(div) for div in divergences)

    for divergence, base_name in zip(divergences, filenames):
        fig, ax = plt.subplots(figsize=(12, 6))
        x = range(0, len(divergence) * step, step)

        # Set color and label based on file name
        if "reverse" in base_name.lower():
            color = 'red'
            label = 'Reverse Divergence'
        elif "forward" in base_name.lower():
            color = 'black'
            label = 'Forward Divergence'
        else:
            color = 'green'
            label = 'Divergence'

        # Plot divergence
        ax.plot(x, divergence, color=color, label=label)

        ax.set_xlabel('Position (bp)', fontname='Arial')
        ax.set_ylabel('Average pairwise divergence', fontname='Arial')
        ax.set_title(base_name, fontname='Arial')

        # Set consistent Y-axis limits
        ax.set_ylim(global_min, global_max)

        # Separate the base from the X-axis
        ax.spines['bottom'].set_position(('outward', 10))

        # Add legend
        ax.legend()

        # Save plots
        plt.savefig(f'Output/{base_name}_divergence.pdf', format='pdf', dpi=300, bbox_inches='tight')
        plt.savefig(f'Output/{base_name}_divergence.png', format='png', dpi=300, bbox_inches='tight')
        plt.close()

def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(description='Compare sequences using dotplot and window similarity')
    parser.add_argument('-dotplot', action='store_true', help='Generate dotplot only')
    parser.add_argument('-window', action='store_true', help='Generate window similarity plot only')
    parser.add_argument('-both', action='store_true', help='Generate both plots')
    
    # Add arguments for file/directory paths
    parser.add_argument('-f', '--fasta', type=str, help='FASTA file for dotplot analysis (must contain exactly 2 sequences)')
    parser.add_argument('-d', '--directory', type=str, help='Directory containing alignment files for window analysis')

    args = parser.parse_args()

    # Create output directory
    setup_output_directory()

    # Perform requested analyses
    if args.dotplot or args.both:
        if not args.fasta:
            print("Error: Please specify a FASTA file using -f or --fasta for dotplot analysis")
            return
        
        if not os.path.exists(args.fasta):
            print(f"Error: FASTA file '{args.fasta}' not found")
            return
        
        try:
            seq1, seq2 = read_sequences(args.fasta)
            matches_forward, matches_reverse = create_dotplot(seq1, seq2)
            plot_dotplot(seq1, seq2, matches_forward, matches_reverse)
            print("Dotplot analysis complete")
        except Exception as e:
            print(f"Error in dotplot analysis: {e}")

    if args.window or args.both:
        if not args.directory:
            print("Error: Please specify a directory using -d or --directory for window analysis")
            return
        
        if not os.path.exists(args.directory):
            print(f"Error: Directory '{args.directory}' not found")
            return
        
        # Process all alignment files in the specified directory
        alignment_files = glob.glob(os.path.join(args.directory, '*.fasta'))
        if not alignment_files:
            print(f"No FASTA files found in directory '{args.directory}'")
        else:
            try:
                plot_window_divergence(alignment_files)
                print("Window divergence analysis complete")
            except Exception as e:
                print(f"Error in window analysis: {e}")

    if not (args.dotplot or args.window or args.both):
        print("Please specify analysis type: -dotplot, -window, or -both")
        print("Use -h or --help for more information")

if __name__ == "__main__":
    main()