# coding=utf-8

from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics.GenomeDiagram._Colors import ColorTranslator
# from Bio.Graphics import GenomeDiagram
from Bio import SeqIO
from reportlab.lib.units import cm
from reportlab.lib import colors
import random
import os
import read_table

codec = {'nt': 'gbk'}.get(os.name, 'utf-8')
abs_dir = os.path.abspath(os.path.join(__file__, os.path.pardir))


def get_info(path):
    info = read_table.read_excel(path)
    info = zip(*info)
    info = [(l[0], [int(i) for i in l[1:] if i != '']) for l in info]
    info_d = {}

    def _add(k, l):
        def _add_v(v):
            if not info_d.has_key(v):
                info_d[v] = {}
            vd = info_d[v]
            if not vd.has_key(k):
                vd[k] = 0
            vd[k] += 1

        [_add_v(v) for v in l]

    categories = zip(*info)[0]
    [_add(k, l) for k, l in info]

    info_d = sorted(info_d.items(), key=lambda x: x[0])
    return info_d, categories


def read_colors():
    with open(os.path.join(abs_dir, 'colors.pref')) as fp:
        lines = [l.strip() for l in fp]
    lines = [l for l in lines if not l.startswith('#') and l != '']
    inner_start = lines.index('[inner]')
    outer_start = lines.index('[outer]')
    if inner_start > outer_start:
        inner = lines[inner_start + 1:]
        outer = lines[outer_start + 1:inner_start]
    else:
        outer = lines[outer_start + 1:]
        inner = lines[inner_start + 1:outer_start]

    def _conv(color):
        if len(color) == 6:
            try:
                color.decode('hex')
                return '#' + color
            except:
                pass
        return color

    inner, outer = [map(_conv, c) for c in (inner, outer)]
    return inner, outer


def draw_by_bio(info, cates, long_intervals, seqlen, filename):
    import diagram as _dia
    reload(_dia)
    diagram = _dia.Diagram('Test Diagram')

    inner_colors, outer_colors = read_colors()
    if len(cates) > len(outer_colors):
        raise Exception('outer_colors num: %s, cannot satisfy category num: %s' % (len(outer_colors), len(cates)))
    color_map = zip(cates, outer_colors)
    color_trans = ColorTranslator()
    color_name_pairs = [(color_trans.translate(t[1]), t[0]) for t in color_map]
    color_map = dict(color_map)

    def _add_feats(feature_set, i, width=3, minimal_angle_margin=20):
        loc, concur_cates = info[i]
        concur_cates = concur_cates.keys()
        concur_num = len(concur_cates)

        # -------- calc width --------
        prev = 0
        if i > 0:
            prev = info[i - 1][0]

        next = seqlen
        if i < len(info) - 1:
            next = info[i + 1][0]

        width = min((next - prev) / concur_num, width)
        if width < 1:
            raise Exception('too narrow')

        # -------- add feature --------
        if concur_num > 1:
            angle_sum = (concur_num - 1) * minimal_angle_margin
            assert angle_sum < 180
            angle_slice = angle_sum / (concur_num - 1)

        def _add_feat(cate_i):
            color = color_map[concur_cates[cate_i]]
            start = loc  # - (cate_i - concur_num / 2) * width
            end = start + width
            feat = SeqFeature(FeatureLocation(start, end, strand=1))
            titled_angle = 0
            if concur_num > 1:
                titled_angle = cate_i * angle_slice - angle_sum / 2
            feature_set.add_feature(feat, color=color, tilted_angle=titled_angle)

        [_add_feat(cate_i) for cate_i in xrange(concur_num)]

    def _add_intervals(feature_set, i):
        start, end = long_intervals[i]
        feat = SeqFeature(FeatureLocation(start, end, strand=-1))
        color = inner_colors[i % len(inner_colors)]
        feature_set.add_feature(feat, color=color)

    def _add(track):
        feature_set = track.new_set()
        [_add_feats(feature_set, i, 1) for i in xrange(len(info))]
        if long_intervals is not None:
            [_add_intervals(feature_set, i) for i in xrange(len(long_intervals))]

    def _track(track_level):
        track = diagram.new_track(track_level, greytrack=False)
        _add(track)

    _track(1)

    diagram.draw(format='circular', circular=True, pagesize=(8 * cm, 8 * cm), fragments=1, orientation='portrait',
                 start=0, end=seqlen, circle_core=.7, color_name_pairs=color_name_pairs)
    diagram.write("%s.svg" % filename, "svg")
    print u'输出为%s.svg' % filename


def single_proc(args):
    args = args.split('  ')
    name, seqlen = args[:2]
    seqlen = int(seqlen)

    long_intervals = None
    if len(args) > 2:
        long_intervals = args[2].split(',')
        long_intervals = [map(int, tup.split('-')) for tup in long_intervals]
    info, cates = get_info(name)
    dst_name = os.path.splitext(name)[0].decode(codec)
    draw_by_bio(info, cates, long_intervals, seqlen, dst_name)


def main():
    while True:
        try:
            opt = raw_input('\n\n1. single file processing\n2. batch script processing\noption (1/2): ')
            if opt == '1':
                args = raw_input(u'输入文件名、序列总长、长区间(start1-end1,start2-end2)，两个空格分开：\n'.encode(codec))
                single_proc(args)
            elif opt == '2':
                path = raw_input(u'输入批处理脚本的路径，批处理脚本的每行与单文件处理格式一致：\n'.encode(codec))
                batch_proc(path)
        except Exception as e:
            if isinstance(e, KeyboardInterrupt):
                raise e
            else:
                print e


def batch_proc(path):
    with open(path) as fp:
        lines = fp.read().splitlines()
    lines = [l.decode('utf-8').encode(codec) for l in lines if l != '']
    [single_proc(l) for l in lines]


def stat():
    data = '''
        CGCTA
        CGCTA
        ATCAC
        ATCAC
        TCCCA
        TAAAG
        GAACC
        GGGCG
        GCAAA
        GCAAA
    '''.splitlines()
    data = [seq.strip() for seq in data]
    map = {'A': 'T',
           'T': 'A',
           'C': 'G',
           'G': 'C'}

    result = {}

    # -------- add --------
    def _add(seq):
        seq = seq.upper()
        result[seq] = result.get(seq, 0) + 1

    [_add(seq) for seq in data if seq != '']

    # -------- detect --------
    def _detect(seq):
        rev_seq = ''.join([map[c] for c in seq])
        if not result.has_key(rev_seq):
            del result[seq]
        else:
            print seq, rev_seq

    [_detect(seq) for seq in result.keys()]
    import ipdb
    ipdb.set_trace()


if __name__ == '__main__':
    main()
