import logging
import os
import pickle
from Bio import SeqIO
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import FeatureLocation, SeqFeature
from reportlab.lib import colors
from .multiplex_scheme import MultiplexScheme

logger = logging.getLogger('Primal Log')

class SchemeReporter(MultiplexScheme):
    """Reporting methods to extend MultiplexScheme"""

    def write_bed(self, path='./'):
        logger.info('Writing BED')
        filepath = os.path.join(path, '{}.bed'.format(self.prefix))
        with open(filepath, 'w') as bedhandle:
            for r in self.regions:
                print >>bedhandle, '\t'.join(map(
                    str, [self.primary_reference.id, r.top_pair.left.start, r.top_pair.left.end, r.top_pair.left.name, r.pool]))
                print >>bedhandle, '\t'.join(map(str, [self.primary_reference.id, r.top_pair.right.end,
                                                       r.top_pair.right.start, r.top_pair.right.name, r.pool]))

    def write_tsv(self, path='./'):
        logger.info('Writing TSV')
        filepath = os.path.join(path, '{}.tsv'.format(self.prefix))
        with open(filepath, 'w') as tsvhandle:
            print >>tsvhandle, '\t'.join(
                ['name', 'seq', 'length', '%gc', 'tm (use 65)'])
            for r in self.regions:
                left = r.top_pair.left
                right = r.top_pair.right
                print >>tsvhandle, '\t'.join(
                    map(str, [left.name, left.seq, left.length, '%.2f' %left.gc, '%.2f' %left.tm]))
                print >>tsvhandle, '\t'.join(
                    map(str, [right.name, right.seq, right.length, '%.2f' %right.gc, '%.2f' %right.tm]))
                if r.alternates:
                    for alt in r.alternates:
                        print >>tsvhandle, '\t'.join(map(str, [alt.name, alt.seq, alt.length, '%.2f' %alt.gc, '%.2f' %alt.tm]))

    def write_pickle(self, path='./'):
        logger.info('Writing pickles')
        filepath = os.path.join(path, '{}.pickle'.format(self.prefix))
        with open(filepath, 'wb') as pickleobj:
            pickle.dump(self.regions, pickleobj)

    def write_refs(self, path='./'):
        logger.info('Writing references')
        filepath = os.path.join(path, '{}.fasta'.format(self.prefix))
        with open(filepath, 'w') as refhandle:
            SeqIO.write(self.references, filepath, 'fasta')

    def write_schemadelica_plot(self, path='./'):
        logger.info('Writing plot')
        gd_diagram = GenomeDiagram.Diagram("Primer Scheme", track_size=1)
        scale_track = GenomeDiagram.Track(
            name='scale', scale=True, scale_fontsize=10, scale_largetick_interval=1000, height=0.1)
        gd_diagram.add_track(scale_track, 2)

        primer_feature_set_1 = GenomeDiagram.FeatureSet()
        primer_feature_set_2 = GenomeDiagram.FeatureSet()

        for r in self.regions:
            cols1 = [self.primary_reference.id, r.top_pair.left.start,
                     r.top_pair.left.end, r.top_pair.left.name, r.pool]
            cols2 = [self.primary_reference.id, r.top_pair.right.end,
                     r.top_pair.right.start, r.top_pair.right.name, r.pool]
            region = str(r.region_num)
            fwd_feature = SeqFeature(FeatureLocation(
                int(cols1[1]), int(cols1[2]), strand=0))
            rev_feature = SeqFeature(FeatureLocation(
                int(cols2[1]), int(cols2[2]), strand=0))
            region_feature = SeqFeature(FeatureLocation(
                int(cols1[1]), int(cols2[2]), strand=0))
            if int(region) % 2 == 0:
                primer_feature_set_1.add_feature(region_feature, color=colors.palevioletred,
                                                 name=region, label=True, label_size=10, label_position="middle", label_angle=0)
                primer_feature_set_1.add_feature(
                    fwd_feature, color=colors.red, name=region, label=False)
                primer_feature_set_1.add_feature(
                    rev_feature, color=colors.red, name=region, label=False)
            else:
                primer_feature_set_2.add_feature(region_feature, color=colors.palevioletred,
                                                 name=region, label=True, label_size=10, label_position="middle", label_angle=0)
                primer_feature_set_2.add_feature(
                    fwd_feature, color=colors.red, name=region, label=False)
                primer_feature_set_2.add_feature(
                    rev_feature, color=colors.red, name=region, label=False)

        primer_track = GenomeDiagram.Track(name="Annotated Features", height=0.1)
        primer_track.add_set(primer_feature_set_1)
        gd_diagram.add_track(primer_track, 4)

        primer_track = GenomeDiagram.Track(name="Annotated Features", height=0.1)
        primer_track.add_set(primer_feature_set_2)
        gd_diagram.add_track(primer_track, 6)

        rows = max(2, int(round(len(self.primary_reference) / 10000.0)))
        gd_diagram.draw(format='linear', pagesize=(300 * rows, 200 * rows),
                        fragments=rows, start=0, end=len(self.primary_reference))

        pdf_filepath = os.path.join(path, '{}.pdf'.format(self.prefix))
        svg_filepath = os.path.join(path, '{}.svg'.format(self.prefix))
        gd_diagram.write(pdf_filepath, 'PDF', dpi=300)
        gd_diagram.write(svg_filepath, 'SVG', dpi=300)
