from reportlab.lib import colors
from reportlab.lib.units import cm
from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
record = SeqIO.read("ar_psb1c3.gb", "genbank")
gd_diagram = GenomeDiagram.Diagram("Yersinia pestis biovar Microtus plasmid pPCP1")
gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
gd_feature_set = gd_track_for_features.new_set()
idx = 0
for feature in record.features:
	if feature.type == "RBS":
		color = colors.red
		feature.strand = None
	elif feature.type == "promoter":
		color = colors.yellow
		feature.strand = None
		idx = feature.location.start	
	elif feature.type == "CDS":
		color = colors.blue
		feature.strand = None
		gd_feature_set.add_feature(feature, sigil="ARROW", 
					   arrowshaft_height=1.0, 
					   label=True)
		continue
	else: 
		continue
	gd_feature_set.add_feature(feature, color=color, label=True)
gd_diagram.draw(format="linear", orientation="landscape", pagesize=(30*cm, 6*cm),
                fragments=1, start=idx, end=len(record))
gd_diagram.write("plasmid_linear.pdf", "PDF")

