#!/usr/bin/env python

"""
HTML interactive sashimi plots 
modified from ggsashimi by guigolab, https://github.com/guigolab/ggsashimi
"""
import re, sys
from argparse import ArgumentParser

import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import pysam


def parse_arguments():
        argparser = ArgumentParser(description="browsable html sashimi")
        argparser.add_argument("-c","--coordinates", type=str, required=True)
        argparser.add_argument("-b", "--bam", type=str, required=True)
        argparser.add_argument("-g", "--gtf", type=str, required=True)
        argparser.add_argument("-o","--output", default="output", type=str)
        argparser.add_argument("-s", "--strand", default="NONE", type=str)
        argparser.add_argument("-v", "--vcf", default=False, type=str)

        return argparser.parse_args()

def parse_coordinates(c):
        c = c.replace(",", "")
        chr = c.split(":")[0]
        start, end = c.split(":")[1].split("-")
        # Convert to 0-based
        start, end = int(start) - 1, int(end)
        return chr, start, end


def count_operator(CIGAR_op, CIGAR_len, pos, start, end, a, junctions):

        # Match
        if CIGAR_op == "M":
                for i in range(pos, pos + CIGAR_len):
                        if i < start or i >= end:
                                continue
                        ind = i - start
                        a[ind] += 1

        # Insertion or Soft-clip
        if CIGAR_op == "I" or CIGAR_op == "S":
                return pos

        # Deletion
        if CIGAR_op == "D":
                pass

        # Junction
        if CIGAR_op == "N":
                don = pos
                acc = pos + CIGAR_len
                if don > start and acc < end:
                        junctions[(don,acc)] = junctions.setdefault((don,acc), 0) + 1

        pos = pos + CIGAR_len

        return pos

def flip_read(s, samflag):
    if s == "NONE" or s == "SENSE":
        return 0
    if s == "ANTISENSE":
        return 1
    if s == "MATE1_SENSE":
        if int(samflag) & 64:
            return 0
        if int(samflag) & 128:
            return 1
    if s == "MATE2_SENSE":
        if int(samflag) & 64:
            return 1
        if int(samflag) & 128:
            return 0


def read_bam(file, chr, start, end , strand):

    # Initialize coverage array and junction dict
    a = {"+" : [0] * (end - start)}
    junctions = {"+": dict()}
    if strand != "NONE":
            a["-"] = [0] * (end - start)
            junctions["-"] = dict()

    samfile = pysam.AlignmentFile(file)

    for read in samfile.fetch(chr, start, end):

        # Move forward if read is unmapped
        if read.is_unmapped:
            continue

        samflag, read_start, CIGAR = read.flag, read.reference_start+1, read.cigarstring

        # Ignore reads with more exotic CIGAR operators
        if any(map(lambda x: x in CIGAR, ["H", "P", "X", "="])):
            continue

        read_strand = ["+", "-"][flip_read(strand, samflag) ^ bool(int(samflag) & 16)]
        
        if strand == "NONE": 
            read_strand = "+"
        CIGAR_lens = re.split("[MIDNS]", CIGAR)[:-1]
        CIGAR_ops = re.split("[0-9]+", CIGAR)[1:]

        pos = read_start

        for n, CIGAR_op in enumerate(CIGAR_ops):
                CIGAR_len = int(CIGAR_lens[n])
                pos = count_operator(CIGAR_op, CIGAR_len, pos, start, end, a[read_strand], junctions[read_strand])

    samfile.close()
    return a, junctions

def parse_gtf(gtf_file, chromosome, start, end):
    annotations = []
    
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            chrom, feature_type, start_pos, end_pos, strand = fields[0], fields[2], int(fields[3])-1, int(fields[4]), fields[6]
            
            if chrom != chromosome or start_pos > end or end_pos < start:
                continue
            
            if feature_type in ["exon"]:
                transcript_id = re.search('transcript_id "([^"]+)"', fields[8]).group(1)
                annotations.append({
                    "chromosome": chrom,
                    "start": start_pos,
                    "end": end_pos,
                    "type": feature_type,
                    "strand": strand,
                    "transcript_id": transcript_id,
                    "attributes": fields[8]
                })
    
    return pd.DataFrame(annotations)

def parse_vcf(vcf_file, chromosome, start, end):
    variants = []
    # Open the VCF file with pysam (handles both gzipped and uncompressed files)
    vcf_in = pysam.VariantFile(vcf_file)

    for record in vcf_in.fetch(chromosome, start, end):
        pos = record.pos
        ref = record.ref
        alts = record.alts  # This is a tuple of alternative alleles

        for alt in alts:
            variants.append((pos, ref, alt))

    return variants

def get_variant_color(alt):
    match alt:
        case 'A':
            return '#43A5BE'
        case 'G':
            return '#4FB06D'
        case 'T':
            return '#F5C26B'
        case 'C':
            return '#F07857'


def create_sashimi(coverage_data, junctions, start, end, annotations, variants):
    # Flatten coverage data into a single array
    positions = list(range(start, end))
    counts = coverage_data["+"]

    df = pd.DataFrame({"Position": positions, "ReadCount": counts})
    

    if variants:
         row_count = 3
         row_heights = [0.6, 0.1, 0.3]
    else:
         row_count = 2
         row_heights = [0.6, 0.4]


    fig = make_subplots(rows=row_count, cols=1, shared_xaxes=True, 
                        row_heights=row_heights,
                        vertical_spacing=0)
    # Create histogram
    bins = end-start
    hist = px.histogram(df, x="Position", y="ReadCount", nbins=bins, color_discrete_sequence=['grey'])
    hist.update_traces(hovertemplate='',hoverinfo='none')
    fig.add_trace(hist.data[0],row=1, col=1)


    colorlist = px.colors.cyclical.Phase

    c = 0
    sign = 1
    max_junction_height = 0
    for (donor, acceptor), count in junctions["+"].items():
        height = np.log(count + 1) * 40
        mid = (acceptor+donor)/2
        text_height = height+2
        color = colorlist[c]
        # Create points for Bezier curve
        bezier_x = np.linspace(donor, acceptor, 99)
        bezier_y = height * 4 * (bezier_x - donor) * (acceptor - bezier_x) / ((acceptor - donor) ** 2)

        fig.add_trace(go.Scatter(x=bezier_x, y=sign * bezier_y, mode='lines', 
                                 line=dict(color=color, width=max(1,int(np.log(count)))), showlegend=False, hovertemplate=f'pos: {donor}-{acceptor}, count: {count}', name=""),
                      row=1, col=1)
        fig.add_annotation(x=mid, y=sign * text_height, showarrow=False, text=count, font=dict(color=color,size=12), bgcolor="white", row=1, col=1)
        max_junction_height = max(max_junction_height,height)
        c+=1
        sign *= -1

    # Group annotations by transcript_id
    grouped_annotations = annotations.groupby("transcript_id")
    y_offset = -5 # Initial y offset for the annotation tracks
    
    for transcript_id, group in grouped_annotations:
        for _, row in group.iterrows():
            if row["strand"] == '-':
                textinshape = "<<<"
            else:
                textinshape = ">>>"
            fig.add_shape(
                type="rect",
                x0=row["start"],
                x1=row["end"],
                y0=y_offset,  # Position each transcript at a different level
                y1=y_offset + 1,
                line=dict(color="Black"),
                fillcolor='Black',
                label=dict(text=textinshape, font=dict(family="Courier New, monospace", size=10, color="White")),
                row=row_count, col=1
            )
            fig.add_trace(go.Scatter(x=[(row["start"] + row["end"]) / 2],
                            y=[y_offset],
                            text=[transcript_id],
                            mode="text", hoverinfo='text',
                            showlegend=False, textfont=dict(color='rgba(0,0,0,0)')),
                        row=row_count, col=1)
        
        y_offset -= 3  # Move to the next level for the next transcript

    if variants:
        variant_positions = [v[0] for v in variants]
        variant_labels = [f"{v[1]}>{v[2]}" for v in variants]
        variant_colors = [get_variant_color(v[2]) for v in variants]

        fig.add_trace(go.Scatter(x=variant_positions,
                        y=[0.5] * len(variant_positions),
                        mode='markers+text',
                        marker_symbol='square',
                        text=variant_labels,
                        hoverinfo='x+text',
                        textposition="top center",
                        hovertemplate='pos: %{x:.0f}<br>variant: %{text}<extra></extra>',
                        showlegend=False,
                        textfont=dict(color='rgba(0,0,0,0)'),
                        marker=dict(color=variant_colors, size=8)),
                    row=2, col=1)

    # Update layout to make space for the annotations
    fig.update_layout(dict1=dict(template="plotly_white"))
    #fig.update_layout(height=800, width=1200)

    hist_y_ticks = np.linspace(-max_junction_height,df["ReadCount"].max(),num=5, dtype=int)
    hist_y_ticks = [int(np.ceil(num / 100.0)) * 100 for num in hist_y_ticks]
    fig.update_yaxes({'tickvals':hist_y_ticks,"ticktext":[t if t>=0 else '' for t in hist_y_ticks]},row=1,col=1)
    fig.update_yaxes(fixedrange=True)
    fig.update_yaxes(showticklabels=False, row=2, col=1)

    
    fig.update_xaxes(
        tickformat=".0f",
        ticksuffix="", 
    )

    fig.update_layout(
        title={
            'text': f'{chr}:{start}-{end}',
            'x': 0.5,
            'xanchor': 'center', 
            'font': {
                'size': 14, 
                'color': 'Grey' 
            }
        },
    )
    if variants:
        fig.update_yaxes(showticklabels=False, row=3, col=1)
        fig.update_xaxes(title_text="Genomic Position", row=3, col=1)
        fig.update_yaxes(title_text="GTF", row=3, col=1)
        fig.update_yaxes(title_text="Variants", row=2, col=1)
    else:
        fig.update_xaxes(title_text="Genomic Position", row=2, col=1)
        fig.update_yaxes(title_text="GTF", row=2, col=1)
    
    fig.update_yaxes(title_text="Counts", row=1, col=1)

    return fig


if __name__ == "__main__":
        
        args = parse_arguments()
        chr, start, end = parse_coordinates(args.coordinates)

        cov, junct = read_bam(args.bam, chr, start, end, args.strand)

        annotations = parse_gtf(args.gtf, chr, start, end)

        if args.vcf:
             variants = parse_vcf(args.vcf, chr, start, end)
        else:
             variants = False

        fig = create_sashimi(cov, junct, start, end, annotations, variants)

        if '.html' in args.output.lower():
             out_name = args.output
        else:
             out_name = f'{args.output}.html'

        fig.write_html(out_name)
    

        
