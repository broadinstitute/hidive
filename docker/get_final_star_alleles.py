'''Determine the best star allele based on site and genomic matches.'''
import argparse
import os
import sys
import csv
import re
import math
from collections import defaultdict, namedtuple
from pathlib import Path
from pysam import *

VcfRec = namedtuple('VcfRec', ['pos', 'ref', 'alt'])
StarMatch = namedtuple('StarMatch', ['star_allele', 'match', 'total'])

def main():
    args = parse_args()
    assert args.pharmvar_dir.is_dir()
    best_genomic_matches = find_best_genomic_matches(
        args.sam, args.exclude_haps, args.ignore_alleles
    )
    sites_per_star_allele = load_star_allele_definitions(
        args.pharmvar_dir, args.seq_start_pos
    )
    all_matched_sites = read_in_all_matched_sites(
        args.matched_sites, args.exclude_haps, args.ignore_alleles
    )

    with open(args.out, 'w') as out_f:
        w = csv.DictWriter(out_f, delimiter='\t', fieldnames=[
            'Hap', 'Final star allele', 'Best genomic match', 
            'Best genomic match %', 'Best site match', 'Best site match %', 
            'Second-best site match (DIFFERENT major star allele)',
            'Second-best site match %',
        ])
        w.writeheader()
        for hap in all_matched_sites:
            write_alleles(
                w, hap, all_matched_sites[hap], sites_per_star_allele,
                [z[0] for z in best_genomic_matches[hap]]
            )


def parse_args():
    '''Parse command line arguments'''
    p = argparse.ArgumentParser(
        description='Identifies the best star allele per haplotype based on '
        'genomic and site matching'
    )
    p.add_argument(
        'sam', type=Path, help='SAM file with star allele reference sequences '
        '(query) aligned to haplotypes (reference). Please use minimap2 '
        'with --eqx and with a high -N.'
    )
    p.add_argument(
        'matched_sites', type=Path, help='File with matched sites, with header '
        'Sample_hap\tStar_allele\tMatched_sites, where matched sites '
        'is a comma-separated list of pos:ref>alt'
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
        '--ignore_alleles', nargs='+',
        help='List of star alleles to ignore, e.g. duplicates '
        '(*1, *3, *9, *27, *34, *106, *115)',
        default=[
            '*1', '*3', '*9', '*27', '*34', '*106', '*115'
        ]
    )
    p.add_argument(
        '--exclude_haps', nargs='+',
        help='List of haplotypes to exclude, e.g. SV haplotypes (none)',
        default=[]
    )
    if len(sys.argv) < 5:
        p.print_help()
        sys.exit(0)
    return p.parse_args()


def find_best_genomic_matches(sam, exclude_haps=[], ignore_alleles=[]):
    '''Find best-matching star allele(s) for each haplotype by scoring CIGAR 
    string over the genomic locus
    '''
    best_genomic_matches = defaultdict(list)
    align_f = AlignmentFile(sam, 'rb') if (
        sam.suffix == '.bam'
    ) else AlignmentFile(sam)
    for r in align_f:
        hap = r.reference_name
        star_allele = r.query_name.replace('CYP2D6', '')
        if hap in exclude_haps or star_allele in ignore_alleles:
            continue
        # genomic match score = num bases match - (mismatch + del + indel)
        score = 0
        for op, op_len in r.cigartuples:
            if op == CEQUAL:
                score += op_len
            elif op in [CDIFF, CINS, CDEL]:
                score -= op_len
        # Reset best matches if star allele is higher-scoring
        if (
            (hap not in best_genomic_matches) or
            (score > best_genomic_matches[hap][0][1])
        ):
            best_genomic_matches[hap] = [(star_allele, score)]
        # Add to best matches if same score
        elif score == best_genomic_matches[hap][0][1] and (
            (star_allele, score) not in best_genomic_matches[hap]
        ):
            best_genomic_matches[hap].append((star_allele, score))
    return best_genomic_matches


def load_star_allele_definitions(vcf_dir, seq_start_pos):
    '''Load star allele definitions from VCF files'''
    sites_per_star_allele = defaultdict(list)
    for fname in vcf_dir.glob('*.vcf'):
        star_allele = fname.stem.replace('CYP2D6_', '*')
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


def read_in_all_matched_sites(
    matched_sites_fname, exclude_haps=[], ignore_alleles=[]
):
    '''Read in matching sites between each star allele and hap'''
    all_matched_sites = defaultdict(lambda: defaultdict(set))
    with open(matched_sites_fname) as f:
        r = csv.DictReader(f, delimiter='\t')
        for line in r:
            hap = line['Sample_hap']
            star_allele = line['Star_allele']
            if hap in exclude_haps or star_allele in ignore_alleles:
                continue
            for site_str in line['Matched_sites'].split(','):
                all_matched_sites[hap][star_allele].add(
                    VcfRec(int(z[0]), z[1], z[2]) for z in
                    re.split(r'[:>]', site_str)
                )
    return all_matched_sites


def write_alleles(
    w, hap, all_matched_sites_hap, sites_per_star_allele, best_genomic_star
):
    '''Write output TSV line combining best genomic match and site match star
    alleles for haplotype.'''
    site_matches = []
    for star_allele in all_matched_sites_hap:
        matched_sites = len(all_matched_sites_hap[star_allele])
        total_sites = len(sites_per_star_allele[star_allele])
        assert matched_sites <= total_sites
        site_matches.append(
            StarMatch(star_allele, matched_sites, total_sites)
        )
    (
        best_site_percentage, best_site_matches,
        second_best_site_percentage, second_best_site_matches
    ) = get_best_site_matches(site_matches)

    final_star_allele = get_final_star_allele(
        best_genomic_star,
        [z.star_allele for z in best_site_matches],
        best_site_percentage
    )
    out_d = {
        'Hap':hap,
        'Final star allele': final_star_allele,
        'Best genomic match': ','.join(best_genomic_star),
        'Best genomic match %': match_frac_str([
            z for z in site_matches if z.star_allele in best_genomic_star
        ]),
        'Best site match': ','.join([
            z.star_allele for z in best_site_matches
        ]),
        'Best site match %': match_frac_str(best_site_matches),
        'Second-best site match (DIFFERENT major star allele)': ','.join([
            z.star_allele for z in second_best_site_matches
        ]),
        'Second-best site match %': match_frac_str(
            second_best_site_matches
        )
    }
    w.writerow(out_d)


def get_final_star_allele(
    best_genomic_stars, best_site_stars, best_site_percentage
):
    final_star_allele = 'Ambiguous'
    if math.isclose(best_site_percentage, 1.0):
        final_star_alleles = set(best_site_stars).intersection(
            best_genomic_stars
        )
        if len(final_star_alleles) == 1:
            final_star_allele = list(final_star_alleles)[0]
    # Has no sites- assign if best genomic match and no adequate site match
    elif '*1.001' in best_genomic_stars:
        final_star_allele = '*1.001'
    return final_star_allele


def get_best_site_matches(sites):
    '''Get star alleles with highest and second-highest % matching and total
    sites for hap, where second-highest is from different major star allele'''
    # decreasing % then # of sites matched
    sites = sorted(
        sites, key=lambda x:(float(x.match) / x.total, x.total), reverse=True
    )
    best_percentage = match_frac(sites[0])
    best_num_sites = sites[0].total
    second_best_percentage = None
    second_best_num_sites = None
    best_site_matches = []
    second_best_site_matches = []  # from different major star allele only
    for site in sites:
        if (
            math.isclose(match_frac(site), best_percentage) and
            site.total == best_num_sites
        ):  # float equals
            best_site_matches.append(site)
        elif (
            site.star_allele.split('.')[0] in
            [z.star_allele.split('.')[0] for z in best_site_matches]
        ):
            continue
        elif second_best_percentage is None:
            second_best_percentage = match_frac(site)
            second_best_site_matches.append(site)
            second_best_num_sites = site.total
        elif (
            math.isclose(match_frac(site), second_best_percentage) and
            site.total == second_best_num_sites
        ):
            second_best_site_matches.append(site)
        else:
            break
    return (
        best_percentage, best_site_matches, second_best_percentage,
        second_best_site_matches
    )


def match_frac(site):
    return float(site.match) / site.total


def match_frac_str(site_list):
    return ','.join([f'{site.match}/{site.total}' for site in site_list])


if __name__ == '__main__':
    main()
