from Bio import SeqIO, SeqFeature
from Bio.SeqFeature import CompoundLocation, FeatureLocation, SeqFeature
import csv
import sys

# Input file specs:
# column 1: feature name; column 2: start location; column 3: end location; column 4: strand
# strand denoted with +1 for plus strand, -1, for minus strand, 0 for None
##### LOCATIONS SHOULD BE *INCLUSIVE*, AS THEY ARE IN GENBANK FILES

def make_gb_feature(name, start_locus, end_locus, strand, feature_type):
    gb_feature_location = FeatureLocation(start_locus, end_locus, strand=strand)
    gb_feature_type = feature_type
    gb_feature = SeqFeature(gb_feature_location, type=gb_feature_type, strand=strand)
    gb_feature.qualifiers['locus_tag'] = name
    return gb_feature

def append_gb_features(base_gb_in, feature_csv, feature_type, outpath):
    base_record = SeqIO.read(open(base_gb_in), 'genbank')
    feature_list = []

    with open(feature_csv, 'r') as feature_list:
        feature_reader = csv.reader(feature_list, delimiter=',')
        next(feature_reader, None) # ignore the header row
        for entry in feature_reader:
            new_feature = make_gb_feature(entry[0], (int(entry[1])-1), int(entry[2]), int(entry[3]), feature_type)
            base_record.features.append(new_feature)

    #Write new IGR features to the genbank
    SeqIO.write(base_record, open(outpath,'w'), 'genbank')

if __name__ == '__main__':
    if len(sys.argv) == 5:
         append_gb_features(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
    else:
         print("Usage: add_gb_features.py base_record.gb feature_table.csv feature_type outpath.gb")
         sys.exit(0)


