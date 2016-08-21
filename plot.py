import sys
from Bio import SeqIO
from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio.SeqFeature import SeqFeature, FeatureLocation

record = list(SeqIO.parse(open(sys.argv[1], 'r'), 'fasta'))[0]

gd_track_for_features = GenomeDiagram.Track(name="Annotated Features", scale=0, scale_largetick_labels=True, scale_largetick_interval=5000, height=1)
gd_diagram = GenomeDiagram.Diagram("Zika Brazil", track_size=0.8)

scale_track = GenomeDiagram.Track(name='scale', scale=1, scale_largetick_labels=True, scale_largetick_interval=5000, height=0.3)
gd_diagram.add_track(scale_track, 1)

primer_feature_set = GenomeDiagram.FeatureSet()

for line in open(sys.argv[2], 'r'):
	cols = line.strip().split()
	strand = (+1 if cols[0].split('_')[-1] == 'L' else -1)
	colour = (colors.red if cols[0].split('_')[-2] == 'in' else colors.blue)
	feature = SeqFeature(FeatureLocation(int(cols[3]), int(cols[4]), strand=strand))
	primer_feature_set.add_feature(feature, color=colour, name=cols[0], label=True, label_size = 4, label_position="start", label_angle = 45)	
	
primer_track = GenomeDiagram.Track(name="Annotated Features")
primer_track.add_set(primer_feature_set)
gd_diagram.add_track(primer_track, 2)
gd_diagram.draw(format='linear', pagesize=(30*cm,8*cm), fragments=1, start=0, end=len(record))
gd_diagram.write(sys.argv[3], "pdf")
