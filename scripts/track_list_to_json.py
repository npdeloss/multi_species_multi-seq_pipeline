import argparse
import seaborn as sns
import json

def main(i, o, f, g, s, l):
    files = None
    if type(i) is list:
        files = i
    else:
        with open(i) as ifile:
            files = [line.strip() for line in ifile.readlines()]
    igv_session = {}
    igv_session['locus'] = l
    genome = {}
    genome['fastaURL'] = f
    genome['indexURL'] = f'{f}.fai'
    annotation = {}
    annotation['name'] = 'Genes'
    annotation['type'] = 'annotation'
    annotation['format'] = 'gtf'
    annotation['nameField'] = 'gene_name'
    annotation['url'] = g
    annotation['indexURL'] = f'{g}.tbi'
    simple_annotation = {}
    simple_annotation['name'] = 'Gene Bodies'
    simple_annotation['type'] = 'annotation'
    simple_annotation['url'] = s
    simple_annotation['searchable'] = True
    igv_session['reference'] = genome
    track_filenames = files
    track_names = [track.split('/')[-1] for track in track_filenames]
    track_colors = sns.hls_palette(len(track_filenames)).as_hex()
    tracks = []
    for name, filename, color in zip(track_names, track_filenames, track_colors):
        track = {}
        track['url'] = filename
        track['name'] = name
        track['color'] = color
        track['autoscaleGroup'] = '1'
        track['height'] = 100
        track['height'] = 100
    igv_session['tracks'] = [simple_annotation, annotation] + tracks
    with open(output.json, 'w') as out:
        out.write(json.dumps(igv_session, sort_keys=True, indent=4))

def parse_arguments():
    parser = argparse.ArgumentParser(description = 'Convert a list of tracks into a JSON file for use with IGV.js')
    parser.add_argument('i', type=str, help='Text file containing a list of tracks')
    parser.add_argument('o', type=str, help='Output JSON file for use with IGV.js')
    parser.add_argument('f', type=str, help='Reference FASTA file. Should be indexed with samtools faidx.')
    parser.add_argument('g', type=str, help='Sorted and bgzipped annotation GTF file. Should be indexed with tabix.')
    parser.add_argument('s', type=str, help='Search BED file. Name field should correspond to name of locus.')
    parser.add_argument('l', type=str, help='Locus that the visualization with default to, in the format chr:start-end')
    args = parser.parse_args()
    return vars(args)

if __name__ == "__main__":
    args = parse_arguments()
    main(**args)
