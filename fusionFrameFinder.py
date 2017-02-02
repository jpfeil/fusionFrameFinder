#!/usr/bin/env python2.7

import argparse
import collections
from cStringIO import StringIO
import csv
import itertools
import logging
import multiprocessing
import re

try:
    import swalign

except ImportError:
    print 'Please install swalign: \n pip install swalign'


def read_fasta(input_file, alphabet):
    """
    This module will accept an input fasta and yield individual sequences
    one at a time.
    """
    regexp_string = ''.join(['[^', alphabet, alphabet.lower(), ']+'])
    nucs_regexp = re.compile(regexp_string)  # Regexp to identify sequence lines
    first_seq = True  # Used to bypass first seq
    seq, seq_id, seq_comments = None, None, None
    for line in input_file:
        line = line.lstrip()
        if len(line) == 0:  # blank line
            continue
        if line[0] == '>':  # id line
            if first_seq:
                first_seq = False
                #  >seq_id comments becomes ['', 'seq_id','comments']
                temp_id = re.split(r'[>\s,]+', line, maxsplit=2)
                seq_id = temp_id[1]
                #  right strip to remove trailing newline character
                if len(temp_id) == 3:
                    seq_comments = temp_id[2].rstrip()
                seq = ''
                continue
            if nucs_regexp.findall(seq):
                seq = re.sub(nucs_regexp, '', seq)
            yield [seq_id, seq_comments, seq.upper()]
            #  >seq_id comments becomes ['', 'seq_id','comments']
            temp_id = re.split(r'[>\s,]+', line, maxsplit=2)
            seq_id = temp_id[1]
            if len(temp_id) == 3:
                seq_comments = temp_id[2].rstrip()
            seq = ''
            continue
        #  Remove whitespaces in sequence string
        line = ''.join(line.split())
        seq += line
    # If seq hasn't been initialized, the input file was empty.
    try:
        seq
    except UnboundLocalError:
        raise RuntimeError('Input peptides file was empty.')
    if seq:
        if nucs_regexp.findall(seq):
            seq = re.sub(nucs_regexp, '', seq)
        yield [seq_id, seq_comments, seq.upper()]


def get_transcriptome_data(infile):
    """
    Parses Gencode transcript FASTA file and returns CDS sequences keyed by the transcript ID

    :param file infile: Gencode FASTA file object
    :return: Dictionary of transcripts keyed by transcript ID and dictionary of transcript IDs keyed by gene ID
    :rtype: tuple
    """
    regex = r"(?P<transcript_id>ENST[0-9A-Z]+.\d+)\|(?P<gene_id>ENSG[0-9A-Z]+.\d+).*CDS:(?P<start>\d+)-(?P<stop>\d+)"
    gene_transcripts = collections.defaultdict(list)
    transcript_cds = {}
    for header, comment, seq in read_fasta(infile, 'ACGT'):
        match = re.search(regex, header)
        if match:
            # Remove the version number on gene ID
            gene_id, version = match.group("gene_id").split('.')
            transcript_id = match.group("transcript_id")
            # GTF is one-based
            start = int(match.group("start")) - 1
            stop = int(match.group("stop"))
            cds = seq[start: stop]

            gene_transcripts[gene_id].append(transcript_id)
            transcript_cds[transcript_id] = cds
    return transcript_cds, gene_transcripts


def read_fusions(fusion_file):
    """
    Reads in gene fusion predictions in modified BEDPE format. In addition to the basic BEDPE features, this
    function requires the fusion junction sequences and HUGO names for the donor and acceptor genes.

    :param file fusion_file: Fusion calls in BEDPE format
    :returns: list of BEDPE namedtuples
    :rtype: list

    Modified BEDPE format
    chrom1:         Chromosome of first feature
    start1:         Zero-based starting position of the first feature or '.'
    end1:           One-based ending position of the first feature -- 5' fusion breakpoint
    chrom2:         Chromosome of second feature
    start2:         Zero-based starting position of the second feature -- 3' fusion breakpoint
    end2:           One-based ending position of thh second feature or '.'
    name:           Hyphenated Ensembl gene IDs (i.e. ENSG00000168172-ENSG00000165731)
    score:          Optional fusion score
    strand1:        Strand of first feature
    strand2:        Strand of second feature
    junctionSeq1:   Fusion junction sequence in first feature
    junctionSeq2:   Fusion junction sequence in second feature
    hugo1:          HUGO name for first feature
    hugo2:          HUGO name for second feature
    """
    BEDPE = collections.namedtuple('BEDPE',
                                   'chrom1, start1, end1, '
                                   'chrom2, start2, end2, '
                                   'name, score, '
                                   'strand1, strand2, '
                                   'junctionSeq1, junctionSeq2, '
                                   'hugo1, hugo2')

    calls = []
    for line in csv.reader(fusion_file, delimiter='\t'):
        if line[0].startswith('#'):
            continue
        try:
            calls.append(BEDPE(*line))
        except TypeError:
            raise ValueError("ERROR: fusion file is malformed.\n{}".format(read_fusions.__doc__))
        break
    return calls


def align_filter(ref, query, mode, fusion_name=''):
    """
    Aligns query to reference CDS sequence using the Smith-Waterman algorithm. Returns None if the
    alignment is clipped at the fusion boundary.

    :param str ref: In-frame reference transcript
    :param str query: Query transcript
    :param str mode: 'donor' or 'acceptor'
    :return: Alignment features
    :rtype: namedtuple
    """
    alignment_stats = collections.namedtuple('AlignStats', 'qstart, qstop, rstart, rstop, insertions, deletions')

    bounds_regex = re.compile(r'Query\s*:\s*(?P<qstart>\d*)\s*\w*\s*(?P<qstop>\d*)\s*[\|\s]*\s*Ref\s*:\s*(?P<rstart>\d*)\s*\w*\s*(?P<rstop>\d*)')
    match_regex = re.compile(r'Matches: \d+\s\((?P<percent>\d*)')

    match = 5
    mismatch = -1
    scoring = swalign.NucleotideScoringMatrix(match, mismatch)
    sw = swalign.LocalAlignment(scoring)
    alignment = sw.align(ref, query)

    # First check that the donor sequence is in frame
    insertions = 0
    deletions = 0
    for chr, num in alignment.cigar:
        if chr == 'I':
            insertions += num

        elif chr == 'D':
            deletions += num

    # Next grab the alignment statistics
    string = StringIO()
    alignment.dump(out=string)
    dump = string.getvalue()
    string.close()

    # If it's not a near perfect match, then the quality of the assembly may not be good
    m = match_regex.search(dump)
    if m:
        percent = int(m.group('percent'))
        if percent < 99:
            # print('Percent matching %d' % percent)
            # print(dump)
            logging.debug('%s: low percent matching %d' % (fusion_name, percent))
            return

    # If the fusion transcript passes these filters, then grab the bounds of the alignment
    s = bounds_regex.search(dump)
    if s:
        qstart = int(s.group('qstart')) - 1    # Make zero-based
        qstop = int(s.group('qstop'))

        # If the end of the fusion transcript doesn't align, then skip this transcript
        if mode == 'donor' and qstop != len(query):
            logging.debug('%s: donor alignment does not include end of sequence' % fusion_name)
            return

        elif mode == 'acceptor' and qstart != 0:
            logging.debug('%s: acceptor alignment does not include start of sequence' % fusion_name)
            # print('Acceptor doesn\'t start at one')
            # print(dump)
            return

        rstart = int(s.group('rstart')) - 1    # Make zero-based
        rstop = int(s.group('rstart'))

        return alignment_stats(qstart, qstop, rstart, rstop, insertions, deletions)

    else:
        return


def scan_frame(reference_start, mode):
    """
    Find the frame of a sequence relative to an in frame sequence

    :param reference_start:
    :return: Number of bases to slice sequence to make it in-frame
    :rtype: int
    """
    # Find the position of the transcript where the predicted donor sequence starts
    in_frame_adjustment = 0
    while (reference_start + in_frame_adjustment) % 3 != 0:
        in_frame_adjustment += 1 if mode is 'donor' else -1

        if abs(in_frame_adjustment) > 3:
            return

    return in_frame_adjustment


def is_in_frame(donor_sequence, donor_reference, acceptor_sequence, acceptor_reference, fusion_name=''):
    """
    Tests whether two fusion sequences combine to form an in-frame fusion transcript

    :param donor_sequence: 5' donor sequence
    :param donor_reference: Gencode CDS transcript
    :param acceptor_sequence: 3' acceptor sequence
    :param acceptor_reference: Gencode CDS transcript
    :return: True if fusion is in frame
    :rtype: Boolean
    """
    donor_align = align_filter(donor_reference, donor_sequence, 'donor', fusion_name=fusion_name)

    if donor_align is None:
        logging.debug('%s: donor did not align' % fusion_name)
        return False

    # If the alignment is out of frame from the reference, then we cannot infer the frame
    if donor_align.insertions % 3 != 0 or donor_align.deletions % 3 != 0:
        logging.debug('%s: donor alignment frameshift' % fusion_name)
        return False

    donor_frame_adjustment = scan_frame(donor_align.rstart, 'donor')
    if donor_frame_adjustment is None:
        return False

    donor_start = donor_align.qstart + donor_frame_adjustment
    in_frame_donor_seq = donor_sequence[donor_start:]

    acceptor_align = align_filter(acceptor_reference, acceptor_sequence, 'acceptor', fusion_name=fusion_name)

    # If the acceptor sequence could not be found, move onto the next transcript
    if acceptor_align is None:
        logging.debug('%s: acceptor did not align' % fusion_name)
        return False

    # If the alignment is out of frame from the reference, then we cannot infer the frame
    if acceptor_align.insertions % 3 != 0 or acceptor_align.deletions % 3 != 0:
        logging.debug('%s: acceptor alignment frameshift' % fusion_name)
        return False

    acceptor_frame_adjustment = scan_frame(acceptor_align.rstop, 'acceptor')
    if acceptor_frame_adjustment is None:
        return False

    acceptor_stop = acceptor_align.qstop + acceptor_frame_adjustment
    in_frame_acceptor_seq = acceptor_sequence[:acceptor_stop]

    if (len(in_frame_donor_seq) + len(in_frame_acceptor_seq)) % 3 == 0:
        return True


def find_frame_star(args):
    """
    This is a workaround for multiprocessing
    """
    return find_frame(*args)


def find_frame(donor_name, donor_sequence, acceptor_name, acceptor_sequence, transcripts, gene_to_transcript_ids):
    fusion_name = '%s-%s' % (donor_name, acceptor_name)
    logging.info('Processing fusion %s' % fusion_name)

    in_frame_fusions = collections.defaultdict(list)
    out_of_frame_fusions = []

    foundAtLeastOne = None
    for donor_transcript_id in gene_to_transcript_ids[donor_name]:
        for acceptor_transcript_id in gene_to_transcript_ids[acceptor_name]:
            try:
                donor_reference = transcripts[donor_transcript_id]
                acceptor_reference = transcripts[acceptor_transcript_id]
                in_frame = is_in_frame(donor_sequence,
                                       donor_reference,
                                       acceptor_sequence,
                                       acceptor_reference,
                                       fusion_name=fusion_name)
                if in_frame is True:
                    foundAtLeastOne = True
                    in_frame_fusions[fusion_name].append('%s-%s' % (donor_transcript_id, acceptor_transcript_id))

            except KeyError:
                continue

    if foundAtLeastOne is None:
        out_of_frame_fusions.append(fusion_name)

    return in_frame_fusions, out_of_frame_fusions


def main(params):
    """
    fusionFrameFinder aligns predicted 5' and 3' fusion sequences to Gencode reference transcripts to find whether
    the frame is preserved. In-frame transcripts are written to an in_frame_<suffix>.txt file and out of frame fusions
    are written to the out_of_frame_<suffix>.txt file.
    """
    logging.basicConfig(level=getattr(logging, params.log_level), format='%(levelname)s: '
                                                                         '%(message)s')

    fusions = read_fusions(params.fusion_file)

    transcripts, gene_to_transcript_ids = get_transcriptome_data(params.transcript_file)

    in_frame_fusions = collections.defaultdict(list)
    out_of_frame_fusions = []

    pool = multiprocessing.Pool(processes=params.CPU)

    # Need to unpack BEDPE object because it can't be pickled
    donor_names = []
    donor_sequences = []
    acceptor_names = []
    acceptor_sequences = []
    for fusion in fusions:
        donor_name, acceptor_name = fusion.name.split('-')
        donor_names.append(donor_name)
        acceptor_names.append(acceptor_name)
        donor_sequences.append(fusion.junctionSeq1)
        acceptor_sequences.append(fusion.junctionSeq2)

    res = pool.map(find_frame_star,
                   itertools.izip(donor_names,
                                  donor_sequences,
                                  acceptor_names,
                                  acceptor_sequences,
                                  itertools.repeat(transcripts),
                                  itertools.repeat(gene_to_transcript_ids)))

    for in_frame, out_frame in res:
        in_frame_fusions.update(in_frame)
        out_of_frame_fusions.extend(out_frame)

    with open('in_frame_%s.txt' % params.suffix, 'w') as f, open('out_of_frame_%s.txt' % params.suffix, 'w') as g:
        f.write('#%s\t%s\n' % ('GeneIDs', 'TranscriptIDs'))
        for fusion_gene, fusion_trancripts in in_frame_fusions.items():
            combined_transcripts = ','.join(fusion_trancripts)
            f.write('%s\t%s\n' % (fusion_gene, combined_transcripts))

        g.write('#GeneIDs\n')
        for fusion_gene in out_of_frame_fusions:
            g.write('%s\n' % fusion_gene)


def run_frame_finder():
    parser = argparse.ArgumentParser(description=main.__doc__)
    parser.add_argument('--transcripts', dest='transcript_file',
                        type=argparse.FileType('r'),
                        help='Path to GENCODE transcript FASTA file',
                        required=True)
    parser.add_argument('--fusions', dest='fusion_file',
                        help='Path to gene fusion file',
                        type=argparse.FileType('r'),
                        required=True)
    parser.add_argument('--suffix', dest='suffix', type=str, default='fusions',
                        help='Suffix for output file names', required=True)
    parser.add_argument('--CPU', dest='CPU', type=int, default=1,
                        help='Number of CPUs to use')
    parser.add_argument('--log_level', dest='log_level', help='The level of logging above which '
                        'messages should be printed.', required=False, choices={'DEBUG', 'INFO',
                                                                                'WARNING', 'ERROR'},
                        default='INFO')
    params = parser.parse_args()

    return main(params)


if __name__ == '__main__':
    run_frame_finder()