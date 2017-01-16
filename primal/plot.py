import sys
from Bio import SeqIO
from reportlab.lib import colors
from reportlab.lib.units import cm
from reportlab.lib.pagesizes import A4, legal, landscape
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation


def main(prefix, reference):
	gd_diagram = GenomeDiagram.Diagram("Zika Brazil", track_size=1)	
	
	"""
	record = SeqIO.read(sys.argv[1], "genbank")
	gd_diagram = GenomeDiagram.Diagram("Zika Brazil", track_size=1)
	gd_track_for_features = gd_diagram.new_track(2, name="Annotated Features", height=1)
	gd_feature_set = gd_track_for_features.new_set()
	for feature in record.features:
		print feature.type
		if feature.type == "mat_peptide":
			gd_feature_set.add_feature(feature, color=colors.blue, label=True, sigil="ARROW", arrowshaft_height=1, label_position="middle", label_angle = 0, label_size = 10, arrowhead_length=0.5)
	
		if 'UTR' in feature.type:
			gd_feature_set.add_feature(feature, color=colors.lightblue, label=True, sigil="ARROW", arrowshaft_height=1, label_position="middle", label_angle = 0, label_size = 10)
	"""

	scale_track = GenomeDiagram.Track(name='scale', scale=True, scale_largetick_interval=1000, height=0.1)
	gd_diagram.add_track(scale_track, 2)

	primer_feature_set_1 = GenomeDiagram.FeatureSet()
	primer_feature_set_2 = GenomeDiagram.FeatureSet()
	with open(sys.argv[2], 'r') as bedfile:
		while True:
			line1 = bedfile.readline()
			if not line1:
				break
			line2 = bedfile.readline()
			cols1 = line1.strip().split()
			cols2 = line2.strip().split()
			region = cols1[3].split('_')[1]
			print cols1
			print cols2
			fwd_feature = SeqFeature(FeatureLocation(int(cols1[1]), int(cols1[2]), strand=0))
			rev_feature = SeqFeature(FeatureLocation(int(cols2[1]), int(cols2[2]), strand=0))
			region_feature = SeqFeature(FeatureLocation(int(cols1[1]), int(cols2[2]), strand=0))
			if int(region) % 2 == 0:
				primer_feature_set_1.add_feature(region_feature, color=colors.palevioletred, name=region, label=True, label_size = 6, label_position="middle", label_angle = 0)
				primer_feature_set_1.add_feature(fwd_feature, color=colors.red, name=region, label=False)
				primer_feature_set_1.add_feature(rev_feature, color=colors.red, name=region, label=False)
			else:
				primer_feature_set_2.add_feature(region_feature, color=colors.palevioletred, name=region, label=True, label_size = 6, label_position="middle", label_angle = 0)
				primer_feature_set_2.add_feature(fwd_feature, color=colors.red, name=region, label=False)
				primer_feature_set_2.add_feature(rev_feature, color=colors.red, name=region, label=False)

		
	primer_track = GenomeDiagram.Track(name="Annotated Features", height=0.1)
	primer_track.add_set(primer_feature_set_1)
	gd_diagram.add_track(primer_track, 4)

	primer_track = GenomeDiagram.Track(name="Annotated Features", height=0.1)
	primer_track.add_set(primer_feature_set_2)
	gd_diagram.add_track(primer_track, 6)


	gd_diagram.draw(format='linear', pagesize='A4', "landscape", fragments=2, start=0, end=len(reference))

	gd_diagram.write(prefix + '.pdf', "pdf")
	gd_diagram.write(prefix + '.png', 'PNG')

if __name__ == "__main__":
	reference = list(SeqIO.parse(open(sys.argv[1], 'r'), 'fasta'))[0]
	main('test', reference)
