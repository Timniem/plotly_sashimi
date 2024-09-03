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

def create_sashimi(coverage_data, junctions, start, end, annotations):
    # Flatten coverage data into a single array
    positions = list(range(start, end))
    counts = coverage_data["+"]

    # Create a DataFrame for plotting
    df = pd.DataFrame({"Position": positions, "ReadCount": counts})
    

    # Create a subplot figure
    fig = make_subplots(rows=2, cols=1, shared_xaxes=True, 
                        row_heights=[0.6, 0.4],
                        vertical_spacing=0.02)
    # Create histogram
    bins = end-start
    hist = px.histogram(df, x="Position", y="ReadCount", nbins=bins, color_discrete_sequence=['grey'])
    hist.update_traces(hovertemplate='',hoverinfo='none')
    fig.add_trace(hist.data[0],row=1, col=1)


    colorlist = px.colors.cyclical.Phase

    c = 0
    for (donor, acceptor), count in junctions["+"].items():
        height = np.log(count + 1) * 40  # Adjust height for better visibility
        width = (acceptor+donor)/2
        text_height = height+15
        # Create points for Bezier curve
        bezier_x = np.linspace(donor, acceptor, 99)
        bezier_y = height * 4 * (bezier_x - donor) * (acceptor - bezier_x) / ((acceptor - donor) ** 2)

        fig.add_trace(go.Scatter(x=bezier_x, y=bezier_y, mode='lines', 
                                 line=dict(color=colorlist[c]), showlegend=False, hoverinfo='none'),
                      row=1, col=1)
        fig.add_trace(go.Scatter(x=[width], y=[text_height], mode='text', text=count, textfont=dict(color=colorlist[c],size=14), showlegend=False, hoverinfo='none'), row=1, col=1)
        c+=1
    
    # Group annotations by transcript_id
    grouped_annotations = annotations.groupby("transcript_id")
    y_offset = -5 # Initial y offset for the annotation tracks
    
    for transcript_id, group in grouped_annotations:
        for _, row in group.iterrows():
            fig.add_shape(
                type="rect",
                x0=row["start"],
                x1=row["end"],
                y0=y_offset,  # Position each transcript at a different level
                y1=y_offset + 1,
                line=dict(color="Black"),
                fillcolor='Black',
                row=2, col=1
            )
            fig.add_trace(go.Scatter(x=[(row["start"] + row["end"]) / 2],
                            y=[y_offset],
                            text=[transcript_id],
                            mode="text", hoverinfo='text',
                            showlegend=False, textfont=dict(color='rgba(0,0,0,0)')),
                        row=2, col=1)
        
        
        y_offset -= 3  # Move to the next level for the next transcript

    # Update layout to make space for the annotations
    fig.update_layout(dict1=dict(template="plotly_white"))
    fig.update_layout(height=800, width=1200)
    fig.update_layout()
    fig.update_yaxes(fixedrange=True)
    fig.update_yaxes(showticklabels=False, row=2, col=1)
    return fig


if __name__ == "__main__":
        
        args = parse_arguments()
        chr, start, end = parse_coordinates(args.coordinates)

        cov, junct = read_bam(args.bam, chr, start, end, args.strand)

        annotations = parse_gtf(args.gtf, chr, start, end)

        fig = create_sashimi(cov, junct,start,end, annotations)

        if '.html' in args.output.lower():
             out_name = args.output
        else:
             out_name = f'{args.output}.html'

        fig.write_html(out_name)

        