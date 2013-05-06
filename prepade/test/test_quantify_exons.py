import unittest
import numpy as np
from Bio.SeqFeature import SeqFeature, FeatureLocation
from prepade.quantifyexons import cigar_to_spans, spans_are_consistent, compare_aln_to_transcript
from itertools import starmap

def spans_to_locations(spans):
    return [ FeatureLocation(*x) for x in spans ]

class QuantifyExonsTest(unittest.TestCase):

    def test_cigar_to_spans(self):
    #                0123456
    # CIGAR_CHARS = 'MIDNSHP=X'

        cases = [
            (([(0, 75)], 3950850), 
             [(3950850, 3950925)]),

            (([(4, 3), (0, 72) ], 3950849), 
             [(3950849, 3950921)]),

            (([ (4, 3), (0, 65), (3, 1512), (0, 7)], 3954408),
             [(3954408, 3954473),
              (3955985, 3955992)]),

            (([(4, 3), (0, 9), (3, 177), (0, 56), (3, 87), (0, 7) ], 11708587),
             [(11708587, 11708596),
              (11708773, 11708829),
              (11708916, 11708923)]),

            (([(0, 66), (2, 2), (0, 34)], 0),
             [(0, 102)])
        ]

        for args, expected in cases:
            got = cigar_to_spans(*args)
            spans = []
            for feature in got:
                spans.append((feature.start,
                              feature.end))

            self.assertEquals(spans, expected)


    def test_spans_are_consistent1(self):
        exon = FeatureLocation(4351909, 4352081)

        spans = spans_to_locations([
            (4352004, 4352081),
            (4352201, 4352224)])

        res = list(spans_are_consistent(exon, spans))
            
        self.assertEquals([True, True], res)

    def test_spans_are_consistent2(self):
        
        # rev: 

        exon = FeatureLocation(121327520, 121327678)

        spans_f = spans_to_locations([
            (121319862, 121319905),
            (121327520, 121327577)])

        spans_r = spans_to_locations([
            (121327576, 121327622),
            (121327624, 121327678)])

        self.assertEquals([True, True], list(spans_are_consistent(exon, spans_f)))
        self.assertEquals([False, False], list(spans_are_consistent(exon, spans_r)))

def transcript_feature(ref, exons):
    exon_locs = spans_to_locations(exons)
    exons = [ SeqFeature(ref=ref, location=FeatureLocation(*x))
              for x in exons ]

    return SeqFeature(sub_features=exons)
        

class TranscriptQuantTest(unittest.TestCase):

    def test_one_exon_one_segment_inside(self):
        self.check([[10, 20]], 
                   [[12, 18]], True, exon_hits=[True])

    def test_one_exon_one_segment_outside(self):
        self.check([[10, 20]],
                   [[8, 22]], True, exon_hits=[True])

    def test_one_exon_one_segment_overlap_left(self):
        self.check([[10, 20]],
                   [[8,  18]], True, exon_hits=[True])

    def test_one_exon_one_segment_overlap_right(self):
        self.check([[10, 20]],
                   [[12, 22]], True, exon_hits=[True])

    def test_one_exon_one_segment_miss_right(self):
        self.check([[10, 20]], [[22, 30]], False, exon_hits=[False])

    def test_one_exon_one_segment_miss_left(self):
        self.check([[10, 20]],
                   [[22, 30]],
                   False, 
                   exon_hits=[False])

    def test_two_exons_two_segments_hit_exact(self):
        self.check([[10, 20], [30, 40]],
                   [[10, 20], [30, 40]],
                   True, exon_hits=[True, True])

    def test_two_exons_two_segments_hit_overlap_left(self):
        self.check([[10, 20], [30, 40]],
                   [[5, 20], [30, 40]],
                   True, exon_hits=[True, True])
        
    def test_two_exons_two_segments_hit_overlap_right(self):
        self.check([[10, 20], [30, 40]],
                   [[10, 20], [30, 45]],
                   True, exon_hits=[True, True])

    def test_two_exons_one_segment_hit_first_exon(self):
        self.check([[10, 20], [30, 40]],
                   [[10, 20]],
                   True, exon_hits=[True, False])

    def test_two_exons_one_segment_hit_first_exon_cross_edge(self):
        self.check([[10, 20], [30, 40]],
                   [[5, 20]],
                   True, exon_hits=[True, False])

    def test_two_exons_one_segment_hit_last_exon(self):
        self.check([[10, 20], [30, 40]],
                   [[30, 40]], True, [False, True])

    def test_two_exons_one_segment_hit_last_exon_cross_edge(self):
        self.check([[10, 20], [30, 40]],
                   [[30, 45]],
                   True, exon_hits=[False, True])

    def test_two_exons_one_segment_miss_left(self):
        self.check([[10, 20], [30, 40]],
                   [[5, 10]],
                   False, exon_hits=[False, False])

    def test_two_exons_one_segment_miss_right(self):
        self.check([[10, 20], [30, 40]],
                   [[40, 50]],
                   False,
                   exon_hits=[False, False])

    def test_two_exons_one_segment_miss_internal(self):
        self.check([[10, 20], [30, 40]], [[22, 28]], False, exon_hits=[False, False],
                   intron_hits=[True])

    def test_two_exons_one_segment_miss_cross_junction_left(self):
        self.check([[10, 20], [30, 40]], [[15, 25]], False, exon_hits=[True, False], intron_hits=[True])

    def test_two_exons_one_segment_miss_cross_junction_right(self):
        self.check([[10, 20], [30, 40]], [[25, 35]], False, exon_hits=[False, True], intron_hits=[True])

    def test_two_exons_one_segment_miss_cross_junction_both(self):
        self.check([[10, 20], [30, 40]], [[15, 35]], False, exon_hits=[True, True], intron_hits=[True])

    def test_two_exons_one_segment_overlaps_whole_gene(self):
        self.check([[10, 20], [30, 40]], [[5, 45]], False, exon_hits=[True, True], intron_hits=[True])

    def test_one_exon_two_segments(self):
        self.check([[10, 30]], [[12, 18], [22, 27]], False, exon_hits=[True], intron_hits=[],
                   introns=np.zeros((0, 2), int),
                   gaps=np.array([[18, 22]], int))


    def check(self, transcript_in, spans_in, decision,
              exon_hits=None,
              intron_hits=None,
              exons=None,
              spans=None,
              gaps=None,
              introns=None):
        m = compare_aln_to_transcript(transcript_in, spans_in)
        if exon_hits is not None:
            np.testing.assert_equal(exon_hits, m.exon_hits)
        if intron_hits is not None:
            np.testing.assert_equal(intron_hits, m.intron_hits)
        if exons is not None:
            np.testing.assert_equal(exons, m.exons)
        if spans is not None:
            np.testing.assert_equal(spans, m.spans)
        if introns is not None:
            np.testing.assert_equal(introns, m.introns)
        if gaps is not None:
            np.testing.assert_equal(gaps, m.gaps)

        self.assertEqual(decision, m.decision)

#2013-04-25 16:48:12 root         INFO         exon is   93463292-93463472
#2013-04-25 16:48:12 root         DEBUG        spans are 93463398-93463473, 93474313-93474338
#2013-04-25 16:48:12 root         DEBUG      overlap is [True, False]
