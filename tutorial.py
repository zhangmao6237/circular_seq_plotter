# coding=utf-8

def top_down():
    from reportlab.lib import colors
    from reportlab.lib.units import cm
    from Bio.Graphics import GenomeDiagram
    from Bio import SeqIO

    record = SeqIO.read("NC_005816.gb", "genbank")

    gd_diagram = GenomeDiagram.Diagram("Yersinia pestis biovar Microtus plasmid pPCP1")
    gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
    gd_feature_set = gd_track_for_features.new_set()

    for feature in record.features:
        if feature.type != "gene":
            # Exclude this feature
            continue
        if len(gd_feature_set) % 2 == 0:
            color = colors.blue
        else:
            color = colors.lightblue
        gd_feature_set.add_feature(feature, color=color, label=True)

    # gd_diagram.draw(format="linear", orientation="landscape", pagesize='A4', fragments=4, start=0, end=len(record))
    # gd_diagram.write("plasmid_linear.pdf", "PDF")

    circle_core = .8
    gd_diagram.draw(format="circular", circular=True, orientation="landscape", pagesize='A4', start=0, end=len(record),
                    circle_core=circle_core)
    gd_diagram.write("plasmid_circular.pdf", "PDF")


def bottom_up():
    from reportlab.lib import colors
    from reportlab.lib.units import cm
    from Bio.Graphics import GenomeDiagram
    from Bio import SeqIO
    record = SeqIO.read("NC_005816.gb", "genbank")

    # Create the feature set and its feature objects,
    gd_feature_set = GenomeDiagram.FeatureSet()
    for feature in record.features:
        if feature.type != "gene":
            # Exclude this feature
            continue
        if len(gd_feature_set) % 2 == 0:
            color = colors.blue
        else:
            color = colors.lightblue
        gd_feature_set.add_feature(feature, color=color, label=True)
        # (this for loop is the same as in the previous example)

    # Create a track, and a diagram
    gd_track_for_features = GenomeDiagram.Track(name="Annotated Features")
    gd_diagram = GenomeDiagram.Diagram("Yersinia pestis biovar Microtus plasmid pPCP1")

    # Now have to glue the bits together...
    gd_track_for_features.add_set(gd_feature_set)
    gd_diagram.add_track(gd_track_for_features, 1)


def minimal_feats():
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio.Graphics import GenomeDiagram
    from reportlab.lib.units import cm

    gdd = GenomeDiagram.Diagram('Test Diagram')
    gdt_features = gdd.new_track(1, greytrack=False)
    gds_features = gdt_features.new_set()

    # Add three features to show the strand options,
    feature = SeqFeature(FeatureLocation(25, 125), strand=+1)
    gds_features.add_feature(feature, name="Forward", label=True)
    feature = SeqFeature(FeatureLocation(150, 250), strand=None)
    gds_features.add_feature(feature, name="Strandless", label=True)
    feature = SeqFeature(FeatureLocation(275, 375), strand=-1)
    gds_features.add_feature(feature, name="Reverse", label=True)

    gdd.draw(format='circular', circular=True, pagesize=(15 * cm, 15 * cm), fragments=1,
             start=0, end=400, circle_core=.8)
    gdd.write("GD_labels_default.pdf", "pdf")


def minimal_short_feats():
    from Bio.SeqFeature import SeqFeature, FeatureLocation
    from Bio.Graphics import GenomeDiagram
    from Bio import SeqIO
    from reportlab.lib.units import cm
    from reportlab.lib import colors
    import random

    # -------- random allocate --------
    def _alloc(bp_num, group_num, group_size=5, offset=0):
        max_num = bp_num / group_size
        idxes = [random.randrange(0, max_num) for _ in xrange(group_num)]
        names = ['test_%s' % i for i in idxes]
        intervals = [(i * group_size + offset, (i + 1) * group_size + offset) for i in idxes]
        return names, intervals

    phase_1, phase_2 = 500, 3000
    names, intervals = _alloc(phase_1, 13, 5, 0), _alloc(phase_2, 10, 5, phase_1)
    names, intervals = zip(names, intervals)
    reduce_func = lambda x, y: list(x) + list(y)
    names, intervals = [reduce(reduce_func, t) for t in (names, intervals)]

    # names = ['t1', 't2']
    # intervals = [(5, 10), (15, 20)]

    record = SeqIO.read("NC_005816.gb", "genbank")

    def draw_by_bio():
        diagram = GenomeDiagram.Diagram('Test Diagram')

        def _track(track_level):
            track = diagram.new_track(track_level, greytrack=False)
            feature_set = track.new_set()
            for name, it in zip(names, intervals):
                feat = SeqFeature(FeatureLocation(*it, strand=1))
                feature_set.add_feature(feat, name=name, label=True, label_angle=90)

            for i, feat in enumerate(record.features[:8]):
                loc = feat.location
                eta = 4.8045
                feat.location = FeatureLocation(int(loc.start / eta), int(loc.end / eta), strand=-1)
                color = colors.blue if i % 2 == 0 else colors.lightblue
                feature_set.add_feature(feat, color=color, label=True)

        _track(1)
        # _track(2)

        diagram.draw(format='circular', circular=True, pagesize=(15 * cm, 15 * cm), fragments=1,
                     start=0, end=phase_1 + phase_2, circle_core=.9)
        diagram.write("GD_labels_shorts_1.pdf", "pdf")

    def draw_by_dfv():
        from dna_features_viewer import GraphicFeature, CircularGraphicRecord
        import matplotlib.pyplot as plt
        _feat = lambda name, it: GraphicFeature(it[0], it[1], +1, name)
        features = [_feat(name, it) for name, it in zip(names, intervals)]
        record = CircularGraphicRecord(phase_1 + phase_2, features)
        import ipdb
        ipdb.set_trace()

    draw_by_bio()


if __name__ == '__main__':
    minimal_short_feats()
