import unittest
from Bio.SeqFeature import SeqFeature, FeatureLocation
from prepade.quantifyexons import cigar_to_spans

class QuantifyExonsTest(unittest.TestCase):

    def test_cigar_to_spans(self):
    #                0123456
    # CIGAR_CHARS = 'MIDNSHP=X'

        cases = [
            (([(0, 75)], 3950850), 
             [(3950850, 3950924)]),

            (([(4, 3), (0, 72) ], 3950849), 
             [(3950849, 3950920)]),

            (([ (4, 3), (0, 65), (3, 1512), (0, 7)], 3954408),
             [(3954408, 3954472),
              (3955985, 3955991)]),

            (([(4, 3), (0, 9), (3, 177), (0, 56), (3, 87), (0, 7) ], 11708587),
             [(11708587, 11708595),
              (11708773, 11708828),
              (11708916, 11708922)])
        ]

        for args, expected in cases:
            got = cigar_to_spans(*args)
            spans = []
            for feature in got.sub_features:
                spans.append((feature.location.start,
                              feature.location.end))

            self.assertEquals(spans, expected)

