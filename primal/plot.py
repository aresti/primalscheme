import os
import sys

from Bio import SeqIO
from reportlab.lib import colors
from reportlab.lib.units import cm
from reportlab.lib.pagesizes import A4, legal, landscape
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation


def plot_schemeadelica(prefix, results, reference, path='./'):
	gd_diagram = GenomeDiagram.Diagram("Primer Scheme", track_size=1)
	scale_track = GenomeDiagram.Track(name='scale', scale=True, scale_fontsize=10, scale_largetick_interval=1000, height=0.1)
	gd_diagram.add_track(scale_track, 2)

	primer_feature_set_1 = GenomeDiagram.FeatureSet()
	primer_feature_set_2 = GenomeDiagram.FeatureSet()

	for r in results:
		cols1 = [reference.id, r.candidate_pairs[0].left.start, r.candidate_pairs[0].left.end, r.candidate_pairs[0].left.name, r.pool]
		cols2 = [reference.id, r.candidate_pairs[0].right.end, r.candidate_pairs[0].right.start, r.candidate_pairs[0].right.name, r.pool]
		region = str(r.region_num)
		print cols1
		print cols2
		fwd_feature = SeqFeature(FeatureLocation(int(cols1[1]), int(cols1[2]), strand=0))
		rev_feature = SeqFeature(FeatureLocation(int(cols2[1]), int(cols2[2]), strand=0))
		region_feature = SeqFeature(FeatureLocation(int(cols1[1]), int(cols2[2]), strand=0))
		if int(region) % 2 == 0:
			primer_feature_set_1.add_feature(region_feature, color=colors.palevioletred, name=region, label=True, label_size=10, label_position="middle", label_angle = 0)
			primer_feature_set_1.add_feature(fwd_feature, color=colors.red, name=region, label=False)
			primer_feature_set_1.add_feature(rev_feature, color=colors.red, name=region, label=False)
		else:
			primer_feature_set_2.add_feature(region_feature, color=colors.palevioletred, name=region, label=True, label_size=10, label_position="middle", label_angle = 0)
			primer_feature_set_2.add_feature(fwd_feature, color=colors.red, name=region, label=False)
			primer_feature_set_2.add_feature(rev_feature, color=colors.red, name=region, label=False)


	primer_track = GenomeDiagram.Track(name="Annotated Features", height=0.1)
	primer_track.add_set(primer_feature_set_1)
	gd_diagram.add_track(primer_track, 4)

	primer_track = GenomeDiagram.Track(name="Annotated Features", height=0.1)
	primer_track.add_set(primer_feature_set_2)
	gd_diagram.add_track(primer_track, 6)

	rows = int(round(len(reference)/5000))
	gd_diagram.draw(format='linear', pagesize=(1000, 250 * rows), fragments=rows, start=0, end=len(reference))

	png_filepath = os.path.join(path, '{}.png'.format(prefix))
	pdf_filepath = os.path.join(path, '{}.pdf'.format(prefix))
	svg_filepath = os.path.join(path, '{}.svg'.format(prefix))
	gd_diagram.write(png_filepath, 'PNG', dpi=300)
	gd_diagram.write(pdf_filepath, 'PDF', dpi=300)
	gd_diagram.write(svg_filepath, 'SVG', dpi=300)
