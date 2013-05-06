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

        a = spans_to_locations([ (12, 18) ])
        m = compare_aln_to_transcript(np.array([[10, 20]], int), a)
        self.assertTrue(m.decision)
        np.testing.assert_equal([[12, 18]], m.spans)
        np.testing.assert_equal(np.zeros((0, 2)), m.gaps)
        np.testing.assert_equal([[10, 20]], m.exons)
        np.testing.assert_equal(np.zeros((0, 2)), m.introns)
        np.testing.assert_equal(0, m.first_exon_hit)
        np.testing.assert_equal(0, m.last_exon_hit)

    def test_one_exon_one_segment_outside(self):

        a = spans_to_locations([ (8, 22) ])
        m = compare_aln_to_transcript(np.array([[10, 20]], int), a)
        self.assertTrue(m.decision)
        np.testing.assert_equal(0, m.first_exon_hit)
        np.testing.assert_equal(0, m.last_exon_hit)


    def test_one_exon_one_segment_overlap_left(self):

        a = spans_to_locations([ (8, 18) ])
        m = compare_aln_to_transcript(np.array([[10, 20]], int), a)
        self.assertTrue(m.decision)
        np.testing.assert_equal(0, m.first_exon_hit)
        np.testing.assert_equal(0, m.last_exon_hit)

    def test_one_exon_one_segment_overlap_right(self):

        a = spans_to_locations([ (12, 22) ])
        m = compare_aln_to_transcript(np.array([[10, 20]], int), a)
        self.assertTrue(m.decision)
        np.testing.assert_equal(0, m.first_exon_hit)
        np.testing.assert_equal(0, m.last_exon_hit)
        

#2013-04-25 16:48:12 root         INFO         exon is   93463292-93463472
#2013-04-25 16:48:12 root         DEBUG        spans are 93463398-93463473, 93474313-93474338
#2013-04-25 16:48:12 root         DEBUG      overlap is [True, False]
