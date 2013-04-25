import unittest
from Bio.SeqFeature import SeqFeature, FeatureLocation
from prepade.quantifyexons import cigar_to_spans, spans_are_consistent

class QuantifyExonsTest(unittest.TestCase):

    def test_cigar_to_spans(self):
    #                0123456
    # CIGAR_CHARS = 'MIDNSHP=X'

        cases = [
            (([(0, 75)], 3950850, 1), 
             [(3950850, 3950924)]),

            (([(4, 3), (0, 72) ], 3950849, 1), 
             [(3950849, 3950920)]),

            (([ (4, 3), (0, 65), (3, 1512), (0, 7)], 3954408, 1),
             [(3954408, 3954472),
              (3955985, 3955991)]),

            (([(4, 3), (0, 9), (3, 177), (0, 56), (3, 87), (0, 7) ], 11708587, 1),
             [(11708587, 11708595),
              (11708773, 11708828),
              (11708916, 11708922)]),

            (([(0, 28), (1, 1), (0, 46)], 15508982, 1),
             [(15508982, 15509055)])
        ]

        for args, expected in cases:
            got = cigar_to_spans(*args)
            spans = []
            for feature in got.sub_features:
                spans.append((feature.location.start,
                              feature.location.end))

            self.assertEquals(spans, expected)


    def test_spans_are_consistent(self):
        exon = SeqFeature(ref='chr1',
                          location=FeatureLocation(4351909, 4352081))

        spans = [ 
            SeqFeature(
                ref='chr1',
                location=FeatureLocation(4352004, 4352081)),
            SeqFeature(
                ref='chr1',
                location=FeatureLocation(4352201, 4352224))]

        res = list(spans_are_consistent(exon, spans))
            
        self.assertEquals([True, True], res)
