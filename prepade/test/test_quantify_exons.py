import unittest
from Bio.SeqFeature import SeqFeature, FeatureLocation
from prepade.quantifyexons import cigar_to_spans, spans_are_consistent
from itertools import starmap

def spans_to_locations(spans):
    return [ FeatureLocation(*x) for x in spans ]

class QuantifyExonsTest(unittest.TestCase):

    def test_cigar_to_spans(self):
    #                0123456
    # CIGAR_CHARS = 'MIDNSHP=X'

        cases = [
            (([(0, 75)], 3950850, 1), 
             [(3950850, 3950925)]),

            (([(4, 3), (0, 72) ], 3950849, 1), 
             [(3950849, 3950921)]),

            (([ (4, 3), (0, 65), (3, 1512), (0, 7)], 3954408, 1),
             [(3954408, 3954473),
              (3955985, 3955992)]),

            (([(4, 3), (0, 9), (3, 177), (0, 56), (3, 87), (0, 7) ], 11708587, 1),
             [(11708587, 11708596),
              (11708773, 11708829),
              (11708916, 11708923)]),

#            (([(0, 28), (1, 1), (0, 46)], 15508982, 1),
#             [(15508982, 15509055)]),

            (([(0, 66), (2, 2), (0, 34)], 0, 1),
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


#2013-04-25 16:48:12 root         INFO         exon is   93463292-93463472
#2013-04-25 16:48:12 root         DEBUG        spans are 93463398-93463473, 93474313-93474338
#2013-04-25 16:48:12 root         DEBUG      overlap is [True, False]
