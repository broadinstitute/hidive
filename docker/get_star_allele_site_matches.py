'''Get star allele site matches via CIGAR string parsing of reference star
allele vs. haplotype alignments'''
import argparse
import sys
import os
from collections import defaultdict, namedtuple
import pysam
from pathlib import Path

VcfRec = namedtuple('VcfRec', ['pos', 'ref', 'alt'])

def main():
    args = parse_args()
    assert args.pharmvar_dir.is_dir()
    sites_per_star_allele = load_star_allele_definitions(
        args.pharmvar_dir, args.seq_start_pos
    )
    matched_sites = defaultdict(lambda: defaultdict(set))

    sam = args.sam

    # Count non-header lines in SAM file
    num_alignments = 0
    with open(sam) as f:
        for line in f:
            if not line.startswith('@'):
                num_alignments += 1
    print(f"Found {num_alignments} alignments")

    align_n = defaultdict(int)
    if num_alignments > 0:
        align_f = pysam.AlignmentFile(sam, 'rb') if (
            sam.suffix == '.bam'
        ) else pysam.AlignmentFile(sam, check_sq=False)
        for i, r in enumerate(align_f):
            if int(r.get_tag('AS')) < args.min_as:
                continue
            k = (r.query_name, r.reference_name)
            align_n[k] += 1
            if i % 10000 == 0:
                print(f"{i}/{num_alignments}")
            process_alignment(
                r, sites_per_star_allele, matched_sites, args.exclude_haps,
                args.indel_buffer, align_n[k]
            )
        write_output(matched_sites, args.out)
    else:
        with open(args.out, 'w') as out_f:
            out_f.write('Sample_hap\tStar_allele\tMatched_sites\n')


def parse_args():
    '''Parse command line arguments'''
    p = argparse.ArgumentParser(
        description='Identifies star allele site matches from reference star '
        'allele vs. haplotype alignments'
    )
    p.add_argument(
        'sam', type=Path, help='SAM file with star allele reference sequences '
        '(query) aligned to haplotypes (reference). Please use minimap2 '
        'with --eqx and with a high -N.'
    )
    p.add_argument(
        'pharmvar_dir', type=Path, help='Directory with PharmVar VCF files'
    )
    p.add_argument('out', type=Path, help='Final star allele TSV')
    p.add_argument(
        '--seq_start_pos', type=int, default=42124499, 
        help='Starting position of star allele sequences on the chromosome '
        '(42124499)'
    )
    p.add_argument(
        '--exclude_haps', nargs='+',
        help='List of haplotypes to exclude, e.g. SV haplotypes (none)',
        default=[]
    )
    p.add_argument(
        '--indel_buffer', type=int,
        help='Number of surrounding bases that must match between star allele '
        'and haplotype to count as a site match for indels (20).',
        default=20
    )
    p.add_argument(
        '--min_as', type=int,
        help='Minimum value for AS tag (DP alignment score) to retain '
        'alignment (9000)', default=9000
    )
    if len(sys.argv) < 4:
        p.print_help()
        sys.exit(0)
    return p.parse_args()


def load_star_allele_definitions(vcf_dir, seq_start_pos):
    '''Load star allele definitions from VCF files.'''
    sites_per_star_allele = defaultdict(list)
    for fname in vcf_dir.glob('*.vcf'):
        star_allele = fname.stem
        offset = 0
        with open(fname) as f:
            for line in f:
                if not line.strip('"').startswith('#'):
                    chrom, pos, _, ref, alt, *_ = line.strip().split('\t')
                    # Position of variant from start of reference star allele
                    pos = int(pos) - seq_start_pos + offset
                    sites_per_star_allele[star_allele].append(
                        VcfRec(pos, ref, alt)
                    )
                    # Account for position shifts due to indels
                    offset += len(alt) - len(ref)
    return sites_per_star_allele


def process_alignment(
    r, sites_per_star_allele, matched_sites, exclude_haps, indel_buffer, align_n
):
    '''Process a single alignment record.'''
    star_allele = r.query_name.replace('*', '_')
    hap = r.reference_name
    if hap in exclude_haps:
        return
    star_allele_sites = sites_per_star_allele[star_allele]
    sites_pos = [z.pos for z in star_allele_sites]
    star_to_hap_pos = {}
    start_hap_pos = None
    for star_pos, hap_pos in r.get_aligned_pairs():  # (ref, hap)
        if start_hap_pos is None and hap_pos is not None:
            start_hap_pos = hap_pos
        if star_pos in sites_pos:
            star_to_hap_pos[star_pos] = hap_pos
    # Sequences for aligned segments
    hap_seq = r.get_reference_sequence().upper()
    star_seq = r.query_sequence.upper()
    for site in star_allele_sites:
        hap_pos = star_to_hap_pos[site.pos]
        if hap_pos is None:
            continue
        # Pos from start of alignment rather than hap seq
        hap_pos -= start_hap_pos
        match_found = check_site_match(
            site, star_seq, hap_pos, hap_seq, indel_buffer, hap, star_allele
        )
        hap_with_n = f"{hap}_{align_n}" if (align_n > 1) else hap
        if match_found:
            matched_sites[hap_with_n][star_allele].add(site)


def check_site_match(
    site, star_seq, hap_pos, hap_seq, buffer, hap, star_allele
):
    '''Check if a site matches between star allele and haplotype.'''
    hap_alt = hap_seq[hap_pos:hap_pos + len(site.alt)]
    if len(site.ref) == len(site.alt):
        return hap_alt == site.alt
    else:  # Handle indels
        surrounding_hap = hap_seq[
            max(0, hap_pos - buffer):hap_pos + len(site.alt) + buffer
        ]
        surrounding_star = star_seq[
            max(0, site.pos - buffer):site.pos + len(site.alt) + buffer
        ]
        return surrounding_hap == surrounding_star


def write_output(matched_sites, out_file):
    '''Write matched sites to the output file.'''
    with open(out_file, 'w') as out_f:
        out_f.write('Sample_hap\tStar_allele\tMatched_sites\n')
        for sample_hap, star_alleles in sorted(matched_sites.items()):
            for star_allele, sites in sorted(star_alleles.items()):
                sites_str = ','.join(
                    f"{s.pos}:{s.ref}>{s.alt}" for s in sorted(sites)
                )
                out_f.write(
                    f"{sample_hap}\t{star_allele.replace('CYP2D6_', '*')}\t"
                    f"{sites_str}\n"
                )

if __name__ == '__main__':
    main()


