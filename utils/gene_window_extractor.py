# conda activate myenv/
# python gene_window_extractor.py gene_table.pkl ref.fa -o output.pkl
# python gene_window_extractor.py gene_table.pkl ref.fa -o output.pkl --source BestRefSeq
# python gene_window_extractor.py gene_table.pkl ref.fa -o output.pkl --min_window_length 50000
import pandas as pd
import pyfaidx
import argparse

def get_intergenic_sequences(df, fasta_file, source=None):
    """Extracts intergenic sequences, grouped by chromosome."""

    if source:
        df = df[df['source'].isin(source)]

    fa = pyfaidx.Fasta(fasta_file)
    results = []

    for strand in ['+', '-']:
        df_strand = df[df['strand'] == strand].sort_values('start')

        for chrom in df_strand['seqid'].unique():  # Iterate through unique chromosomes
            df_chrom_strand = df_strand[df_strand['seqid'] == chrom]

            for i in range(1, len(df_chrom_strand) - 1):
                prev_gene = df_chrom_strand.iloc[i - 1]
                current_gene = df_chrom_strand.iloc[i]
                next_gene = df_chrom_strand.iloc[i + 1]

                seq_start = prev_gene['end'] + 1
                seq_end = next_gene['start'] - 1

                if seq_start > seq_end:
                    continue  # Skip if start is after end

                try:
                    seq = fa[chrom][seq_start - 1:seq_end].seq
                except KeyError:
                    print(f"Warning: Chromosome {chrom} not found in FASTQ file. Skipping.")
                    continue
                except IndexError:
                    print(f"Warning: Indexing error for gene {current_gene['gene_name']} (start: {seq_start}, end: {seq_end}) on chromosome {chrom}. Check gene coordinates and FASTQ file.")
                    continue

                results.append({
                    'name': current_gene['name'],
                    'chrom': chrom,
                    'gene_start': current_gene['start'],
                    'gene_stop': current_gene['end'],
                    'seq_start': seq_start,
                    'seq_stop': seq_end,
                    'strand': strand,
                    'seq': seq,
                    'source': current_gene['source']
                })

    return pd.DataFrame(results)


def get_gene_window_sequences(df, fasta_file, source=None, max_window_length=178000):
    """
    Extracts gene window sequences centered on the TSS.
    The TSS is determined from associated mRNA features:
      - For '+' strand: the minimum start among mRNA records.
      - For '-' strand: the maximum end among mRNA records.
    If no mRNA records exist, the gene's own start (for '+') or end (for '-') is used.
    mRNA features are grouped with their corresponding gene by matching the mRNA's parent_id with the gene's id.
    The window is constructed to have exactly max_window_length bases (if possible).
    """
    if source:
        df = df[df['source'].isin(source)]
    
    # Load the FASTA file.
    fa = pyfaidx.Fasta(fasta_file)
    results = []
    
    # Separate gene and mRNA features.
    gene_features = df[df['type'] == 'gene'].copy()
    mrna_features = df[df['type'] == 'mRNA'].copy()
    
    for _, gene_info in gene_features.iterrows():

        gene_id = gene_info['id']
        chrom = gene_info['seqid']
        strand = gene_info['strand']
        gene_start = gene_info['start']
        gene_end = gene_info['end']
        gene_length = gene_end - gene_start + 1
        
        # Get mRNA features associated with the gene (matching parent_id to gene id)
        gene_mrnas = mrna_features[mrna_features['parent_id'] == gene_id]
        
        # Determine the TSS.
        if not gene_mrnas.empty:
            if strand == '+':
                tss = gene_mrnas['start'].min()
            else:
                tss = gene_mrnas['end'].max()
        else:
            # Fallback: use gene's own boundary.
            tss = gene_start if strand == '+' else gene_end
        
        # Compute initial window boundaries centered on the TSS.
        half_window = max_window_length // 2
        seq_start = tss - half_window
        seq_end = seq_start + max_window_length - 1  # Ensures fixed window length.
        
        # Get chromosome length.
        chrom_length = len(fa[chrom])

        # Calculate how many bases we need to pad on each side.
        left_pad = 0
        right_pad = 0
        
        # If the left edge goes past the beginning of the chromosome:
        if seq_start < 1:
            left_pad = 1 - seq_start
            seq_start = 1

        # If the right edge goes past the end of the chromosome:
        if seq_end > chrom_length:
            right_pad = seq_end - chrom_length
            seq_end = chrom_length
        
        try:
            # Extract the in-bounds segment. Note that fa[chrom] is zero-indexed in Python,
            # but TSS is presumably 1-based, so we subtract 1 for slicing.
            seq_fragment = fa[chrom][seq_start - 1:seq_end].seq

            # Pad with 'N' to ensure the final sequence has length max_window_length.
            seq = ('N' * left_pad) + seq_fragment + ('N' * right_pad)

            # Confirm we have the fixed length.
            assert len(seq) == max_window_length

        except KeyError:
            print(f"Warning: Chromosome {chrom} not found in FASTA file. Skipping.")
            continue
        except IndexError:
            print(f"Warning: Indexing error for gene {gene_id} (start: {seq_start}, end: {seq_end}) on chromosome {chrom}. Check gene coordinates and FASTA file.")
            continue
        
        results.append({
            'gene_id': gene_id,
            'name': gene_info['name'],
            'chrom': chrom,
            'gene_start': gene_start,
            'gene_end': gene_end,
            'gene_length': gene_length,
            'seq_start': seq_start,
            'seq_end': seq_end,
            'seq_length': len(seq),
            'strand': strand,
            'seq': seq,
            'source': gene_info['source'],  
        })
    
    return pd.DataFrame(results)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Extract gene window or intergenic sequences.")
    parser.add_argument("gene_table", help="Path to the gene table (PKL file).")
    parser.add_argument("fasta_file", help="Path to the FASTA file.")
    parser.add_argument("--source", nargs='+', help="Optional list of sources to filter by.")
    parser.add_argument("-o", "--output", required=True, help="Path to output PKL file.")
    parser.add_argument("--min_window_length", type=int, help="Minimum window length for gene window sequences.")

    args = parser.parse_args()

    df = pd.read_pickle(args.gene_table)

    seqs_df = get_intergenic_sequences(df, args.fasta_file, args.source)

    seqs_df.to_pickle(args.output)
    print(f"Sequences saved to {args.output}")