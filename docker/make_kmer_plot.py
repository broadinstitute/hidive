import sys
import argparse
import json
from pathlib import Path
from collections import defaultdict
import pysam
from yattag import Doc

def main():
    args = parse_args()
    with open(args.plot_options) as f:
        plot_options = json.load(f)
    main_genes, colors, kmer_shape, sort_by = (plot_options[z] for z in [
        'main_genes', 'colors', 'kmer_shape', 'sort_by'
    ])
    kmer_windows, alignments = process_sam(args.sam)
    out_lines = get_out_lines(
        alignments, kmer_windows, main_genes, colors, kmer_shape, sort_by
    )
    output_plot(args.out, out_lines, main_genes, colors, kmer_shape)


def parse_args():
    '''Parse command line arguments'''
    p = argparse.ArgumentParser(
        description='Make a figure with sample sequences divided into '
        'k-mers colored according to which gene they map most closely to.'
    )
    p.add_argument(
        'sam', type=Path, help='Path to sam/ bam of k-mers for each '
        'sample mapped to genes/ sequences of interest. k-mer fasta should '
        'be formatted as <sample>_<start_index>_to_<end_index>'
    )
    p.add_argument(
        'plot_options', type=Path, help='JSON containing plotting parameters. '
        'Keys are "main_genes", "colors", "kmer_shape", and "sort_by" (can be '
        '"colored_blocks" or "sample").'
    )
    p.add_argument('out', type=Path, help='HTML plot output path.')
    if len(sys.argv) < 4:
        p.print_help()
        sys.exit(0)
    return p.parse_args()


def process_sam(sam_fname):
    ''' Read in set of all k-bp windows and samples from sam'''
    alignments = defaultdict(lambda: defaultdict(set))
    max_window = 0
    window_size = None
    align_f = pysam.AlignmentFile(sam_fname, 'rb') if (
        sam_fname.suffix == '.bam'
    ) else pysam.AlignmentFile(sam_fname)
    for r in align_f:
        qname_l = r.query_name.split('_')
        sample = '_'.join(qname_l[:-3])
        window_start = int(qname_l[-3])
        window_end = int(qname_l[-1])
        if not window_size:
            window_size = window_end - window_start + 1
        max_window = max(max_window, window_start)
        if window_start not in alignments[sample]:
            alignments[sample][window_start] = set()
        if r.flag & 0x4 == 0: # is mapped
            gene = r.reference_name.split('_')[0]  # assume seq name starts w/ gene
            # Highest mapq = best-mapping gene
            alignments[sample][window_start].add((gene, int(r.mapq)))
    all_windows = range(0, max_window, window_size)
    return all_windows, alignments


def get_out_lines(
    alignments, kmer_windows, main_genes, colors, kmer_shape='▉',
    sort_by='colored_blocks'
):
    '''Make list with output HTML string for each sample'''
    out_lines = []
    for sample in alignments:
        out_str = ''
        for window in kmer_windows:
            color = assign_kmer_color(
                alignments, sample, window, main_genes, colors
            )
            out_str += f'<span style="color: {color};">{kmer_shape}</span>'
        out_lines.append((sample, out_str))
    if sort_by == 'colored_blocks':
        # sort by most mapped blocks for main genes
        out_lines = sorted(out_lines, key=lambda x:sum([
            x[1].count(colors[z]) for z in main_genes + ['Both']
        ]), reverse=True)
    else:
        assert sort_by == 'sample'
        out_lines = sorted(out_lines, key=lambda x: x[0])
    return out_lines


def assign_kmer_color(alignments, sample, window, main_genes, colors):
    '''Assign k-mer a color based on closest-mapping gene'''
    best_mapping_genes = []
    color_key = 'None'
    if window in alignments[sample]:
        alignments = alignments[sample][window]
        if alignments:
            max_gene_score = max([z[1] for z in alignments])
            best_mapping_genes = [
                z[0] for z in alignments if z[1] == max_gene_score
            ]
            if all([z in best_mapping_genes for z in main_genes]):
                color_key = 'Both'
            elif best_mapping_genes:
                # Prefer main gene(s) as best-mapping
                best_mapping_main = set(best_mapping_genes).intersection(main_genes)
                if best_mapping_main:
                    best_mapping_genes = list(best_mapping_main)
                if len(best_mapping_genes) > 1:
                    print(
                        f'Warning: Ambiguous best-mapping gene for '
                        f'{sample}, index {window}: {", ".join(best_mapping_genes)}'
                    )
                color_key = best_mapping_genes[0]
    else:
        color_key = 'Past_end'
    return colors[color_key]


def output_plot(out_fname, out_lines, main_genes, colors, kmer_shape='▉'):
    '''Create HTML output file containing k-mer plot'''
    doc, tag, text = Doc().tagtext()
    doc.asis('<!DOCTYPE html>')
    with tag('html'):
        with tag('head'):
            with tag('title'):
                text(f"{'/'.join(main_genes)} Mappings")
            with tag('style', type='text/css'):
                text('body { font-family: "DejaVu Sans Mono", monospace; }')
                text('th, td { padding-right: 8px; text-align: left; }')
        with tag('body'):
            with tag('table'):
                with tag('thead'):
                    with tag('tr'):
                        with tag('th'):
                            text('Sample')
                        with tag('th'):
                            text('')
                with tag('tbody'):
                    for out_line in out_lines:
                        with tag('tr'):
                            with tag('td', style='white-space: nowrap;'):
                                doc.asis(f'{out_line[0]}')
                            with tag('td'):
                                doc.asis(str(out_line[1]))
                        text('\n')
        doc.stag('hr')
        doc.asis('<b>K-mer Maps Best To:</b>')
        with tag('table'):
            for gene in colors:
                color = colors[gene]
                with tag('tr'):
                    with tag('td'):
                        if gene == 'Both':
                            text(' and '.join(main_genes))
                        else:
                            text(gene)
                    with tag('td', style=f'color: {color};'):
                        text(kmer_shape)
    with open(out_fname, 'w') as out_f:
        out_f.write(doc.getvalue())


if __name__ == '__main__':
    main()
