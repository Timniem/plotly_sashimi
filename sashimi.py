#!/usr/bin/env python

"""
HTML interactive sashimi plots 
modified from ggsashimi by guigolab, https://github.com/guigolab/ggsashimi
"""
import re, os, subprocess
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
        argparser.add_argument("-o","--output", default="sashimi_output", type=str)
        argparser.add_argument("-s", "--strand", default="NONE", type=str)
        #Optional variants track arguments
        argparser.add_argument("-v", "--vcf", default=False, type=str)
        #Optional spliceAI arguments
        argparser.add_argument("-sa", "--spliceai", default=False, type=bool)
        argparser.add_argument("-gb", "--genomebuild", default="grch37",choices=["grch37","grch38"])
        argparser.add_argument("-t", "--temp", default="./temp/", type=str )
        argparser.add_argument("-r", "--reference", type=str)

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
    start -= 10000
    end += 10000
    
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
    splice_ai = []
    # Open the VCF file with pysam (handles both gzipped and uncompressed files)
    vcf_in = pysam.VariantFile(vcf_file)
    
    for record in vcf_in:
        pos = record.pos
        ref = record.ref
        alts = record.alts  # This is a tuple of alternative alleles
        if 'SpliceAI' in record.info:
                scores = record.info['SpliceAI'][0].split('|')
                #tuple with bases from alt-position and delta scores
                scores_dict = {"AG":(int(scores[6]),float(scores[2])),
                               "AL":(int(scores[7]),float(scores[3])),
                               "DG":(int(scores[8]),float(scores[4])),
                               "DL":(int(scores[9]),float(scores[5]))
                                }
                splice_ai.append((record.pos, scores_dict))
                
        for alt in alts:
            variants.append((pos, ref, alt))

    return variants, splice_ai

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

def subset_vcf(tmp_dir, input_vcf, chrom, start, end):
    """Subset the input VCF to a specific genomic range."""

    tmp_path = os.path.join(tmp_dir, "intermediates")
    output_name = f'{chrom}-{start}-{end}_{os.path.basename(input_vcf)}'
    output_vcf = f'{os.path.join(tmp_path, os.path.splitext(output_name)[0])}'

    if not os.path.exists(tmp_path):
        os.makedirs(tmp_path)
    with pysam.VariantFile(input_vcf) as vcf_in, open(output_vcf, "w") as vcf_out:
        # Copy header
        vcf_out.write(str(vcf_in.header))
        # Filter records by region
        for record in vcf_in.fetch(chrom, start, end):
            vcf_out.write(str(record))

    print(f'Subset VCF {chrom}:{start}-{end} written to {output_vcf}')
    return output_vcf

def annotate_with_spliceai(input_vcf, reference_genome, annotation):
    spliceai_output_name = f'{os.path.splitext(input_vcf)[0]}_spliceai.vcf'
    """Annotate variants using SpliceAI."""
    cmd = [
        "spliceai",
        "-I", input_vcf,
        "-O", spliceai_output_name,
        "-R", reference_genome,  # Reference genome file path
        "-A", annotation  # Output format can be vcf or tsv
    ]
    subprocess.run(cmd)
    print(f"SpliceAI annotated VCF written to {spliceai_output_name}")
    return spliceai_output_name


def create_sashimi(coverage_data, junctions, start, end, annotations, variants, spliceai):
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
    hist = px.histogram(df, x="Position", y="ReadCount", nbins=bins, color_discrete_sequence=['Black'])
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
        bezier_y = (height * 4 * (bezier_x - donor) * (acceptor - bezier_x) / ((acceptor - donor) ** 2))
        steeper = 1 - ((bezier_x - mid) * (mid - bezier_x) / (acceptor - donor) ** 2) * 4
        bezier_y =  bezier_y * steeper

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
        if group['strand'].iloc[0] == '-':
             marker = 'y-left'
        else:
             marker = 'y-right'
        min_transcript = np.inf
        max_transcript = 0
        for _, row in group.iterrows():
            min_transcript = min(min_transcript, row['end'])
            max_transcript = max(max_transcript, row['end'])
            if row["end"] < start:
                 continue
            if row["end"] > end:
                 continue
            fig.add_shape(
                type="rect",
                x0=row["start"],
                x1=row["end"],
                y0=y_offset,  # Position each transcript at a different level
                y1=y_offset - 1,
                line=dict(color="Black"),
                fillcolor='Black',
                row=row_count, col=1
            )
            fig.add_trace(go.Scatter(x=[(row["start"] + row["end"]) / 2],
                            y=[y_offset],
                            text=[transcript_id],
                            mode="text", hoverinfo='text',
                            showlegend=False, textfont=dict(color='rgba(0,0,0,0)')),
                        row=row_count, col=1)
        
        if min_transcript < start:
            min_transcript = start
        if max_transcript > end:
             max_transcript = end
        x_line = np.arange(min_transcript, max_transcript,step=100)
        fig.add_trace(go.Scatter(x=x_line, y=[y_offset-.5]*len(x_line),
                                 mode='lines+markers', marker_line_color='black',
                                  showlegend=False, marker_color="black", marker_symbol=marker,
                                  marker_size=5,marker_line_width=0.5, hovertemplate='',hoverinfo='none'), row=row_count,col=1)
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
        fig.add_hline(y=0.5, line_color="grey",opacity=0.2, row=2, col=1)
        
        if spliceai:
            fig.add_hline(y=3, line_color="grey", line_dash="dot", opacity=0.2, row=2, col=1)
            fig.add_hline(y=5, line_color="grey", line_dash="dot", annotation_text="SpliceAI prediction", 
                            annotation_position="bottom right", opacity=0.2, row=2, col=1)
            spliceai_positions = [s[0] for s in spliceai]
            spliceai_scores = [s[1] for s in spliceai]
            for position, score_dict in zip(spliceai_positions, spliceai_scores):
                for metric in score_dict.keys():
                    if score_dict[metric][1] > 0:
                        if "A" in metric:
                             color = 'orange'
                             ypos = 3
                        else:
                             ypos = 5
                             color = 'blue'
                        if "L" in metric:
                             ydir = - score_dict[metric][1] * 2
                        else:
                             ydir = + score_dict[metric][1] * 2
                        fig.add_shape(type='rect',
                                    x0=position + score_dict[metric][0] - .5,
                                    x1=position + score_dict[metric][0] + .5,
                                    y0=ypos,  # Position each transcript at a different level
                                    y1=ypos+ydir,
                                    showlegend=False,
                                    fillcolor=color,
                                    line=dict(color=color, width=4),
                                    row=2, col=1)
                        fig.add_trace(go.Scatter(x=[position + score_dict[metric][0]],
                                        y=[ypos],
                                        text=metric,
                                        mode='markers+text',
                                        textposition="top center",
                                        hovertemplate=f'SpliceAI - pos: {position + score_dict[metric][0]} - {metric}: {score_dict[metric][1]}',
                                        name="",
                                        marker=dict(color='rgba(0,0,0,0)'),
                                        showlegend=False, textfont=dict(color='rgba(0,0,0,0)')),
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
             subset_vcf_name = subset_vcf(args.temp, args.vcf, chr, start, end)
             if args.spliceai:
                spliceai_vcf_name = annotate_with_spliceai(subset_vcf_name, args.reference, args.genomebuild)
                variants, spliceai = parse_vcf(spliceai_vcf_name, chr, start, end)
             else:
                variants, spliceai = parse_vcf(subset_vcf_name, chr, start, end)
        else:
             variants, spliceai = False, False

        fig = create_sashimi(cov, junct, start, end, annotations, variants, spliceai)

        if '.html' in args.output.lower():
             out_name = args.output
        else:
             out_name = f'{args.output}.html'

        fig.write_html(out_name)
    

        
