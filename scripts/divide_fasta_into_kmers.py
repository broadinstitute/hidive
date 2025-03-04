import sys
import argparse
from pathlib import Path

def main():
    args = parse_args()
    first_line = True
    with open(args.out, 'w') as out_f, open(args.fasta) as f:
        for line in f:
            if line.startswith('>'):
                sample_name = line.strip()
                i = 0
            else:
                for ch in line.strip():
                    if i % args.kmer_size == 0:
                        if not first_line:
                            out_f.write('\n')
                        first_line = False
                        out_f.write(
                            f'{sample_name}_{i}_to_{i + args.kmer_size - 1}\n'
                        )
                    out_f.write(ch)
                    i += 1

def parse_args():
    '''Parse command line arguments'''
    p = argparse.ArgumentParser(
        description='Divide input fasta sequences into k-mer windows.'
    )
    p.add_argument('fasta', type=Path, help='Input fasta')
    p.add_argument('out', type=Path, help='Output fasta')
    p.add_argument(
        '-k', '--kmer_size', type=int, help='K-mer size in base pairs (100)',
        default=100
    )
    if len(sys.argv) < 3:
        p.print_help()
        sys.exit(0)
    return p.parse_args()


if __name__ == '__main__':
    main()