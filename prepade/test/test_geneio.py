import unittest
import numpy as np
from prepade.geneio import read_exons, ExonIndex, parse_gtf_to_genes

class ExonIndexTest(unittest.TestCase):
    
    def test_exon_index(self):

        exons = ['chr1:5-10',
                 'chr1:11-20',
                 'chr1:15-30',
                 'chr1:18-25',
                 'chr1:35-40']

        idx = ExonIndex(read_exons(exons))
        
        np.testing.assert_almost_equal(
            idx.start['chr1'],
            [4, 10, 14, 17, 20, 25, 30, 34, 40])
            
        self.assertEquals(
            [set([('chr1:5-10', 'chr1',  4, 10)]),
             set([('chr1:11-20', 'chr1', 10, 20)]),
             set([('chr1:11-20', 'chr1', 10, 20), ('chr1:15-30', 'chr1', 14, 30)]),
             set([('chr1:11-20', 'chr1', 10, 20), ('chr1:15-30', 'chr1', 14, 30), ('chr1:18-25', 'chr1', 17, 25)]),
             set([('chr1:15-30', 'chr1', 14, 30), ('chr1:18-25', 'chr1', 17, 25)]),
             set([('chr1:15-30', 'chr1', 14, 30)]),
             set(),
             set([('chr1:35-40', 'chr1', 34, 40)]),
             set()],            
            idx.keys['chr1'])

        self.assertEquals(0, len(list(idx.get_exons('chr1', 0, 3))))
        self.assertEquals(0, len(list(idx.get_exons('chr1', 0, 4))))
        self.assertEquals(1, len(list(idx.get_exons('chr1', 0, 5))))
        self.assertEquals(1, len(list(idx.get_exons('chr1', 3, 5))))
        self.assertEquals(1, len(list(idx.get_exons('chr1', 3, 7))))
        self.assertEquals(2, len(list(idx.get_exons('chr1', 7, 12))))
        self.assertEquals(3, len(list(idx.get_exons('chr1', 7, 17))))
        self.assertEquals(2, len(list(idx.get_exons('chr1', 10, 17))))
        self.assertEquals(1, len(list(idx.get_exons('chr1', 28, 29))))
        self.assertEquals(0, len(list(idx.get_exons('chr1', 30, 34))))
        self.assertEquals(1, len(list(idx.get_exons('chr1', 35, 45))))
        self.assertEquals(0, len(list(idx.get_exons('chr1', 45, 100))))


    def test_multiple_refs(self):

        exons = [
            'chr1:5-10',
            'chr1:11-20',
            'chr1:15-30',
            'chr1:18-25',
            'chr1:35-40',

            'chr1:5-10',
            'chr2:11-20',
            'chr3:15-30',
            'chr4:18-25',
            'chr3:35-40',

            'chr1:5-10',
            'chr3:11-20',
            'chr13:15-30',
            'chr3:318-25',
            'chr3:35-40',
        ]
        idx = ExonIndex(read_exons(exons))

class GtfParserTest(unittest.TestCase):

    def setUp(self):
        self.parsed = list(parse_gtf_to_genes('prepade/test/mm9.gtf'))

    def test_parsed_transcript_ids(self):
        transcripts = self.parsed
        self.assertEquals([ t.id for t in self.parsed ], [
            'NM_001011874(refseq)::::uc007aeu.1(ucsc)::::OTTMUST00000065166(vega)',
            'NM_001177658(refseq)::::uc007aff.2(ucsc)',
            'NM_001195662(refseq)::::uc007aew.1(ucsc)',
            'NM_011283(refseq)::::uc007aex.2(ucsc)',
            'NM_011441(refseq)::::uc007aez.1(ucsc)',
            'NM_025300(refseq)::::uc007afd.2(ucsc)',
            'NM_026303(refseq)::::uc009oat.1(ucsc)',
            'NR_027988(refseq)::::uc009oaq.2(ucsc)',
            'NR_033530(refseq)::::uc007afe.2(ucsc)',
            'NR_039546(refseq)',
            'OTTMUST00000065165(vega)',
            'OTTMUST00000086624(vega)',
            'OTTMUST00000086625(vega)',
            'uc007aet.1(ucsc)',
            'uc007aev.1(ucsc)',
            'uc007aey.1(ucsc)',
            'uc007afa.1(ucsc)',
            'uc007afb.1(ucsc)',
            'uc007afc.1(ucsc)',
            'uc009oap.1(ucsc)',
            'uc009oar.2(ucsc)',
            'uc009oas.1(ucsc)'])

    def test_parsed_refs(self):
        self.assertEquals(set(['chr1', 'chr9']), 
                          set([ t.ref for t in self.parsed ]))

    def test_parsed_starts(self):
        got_starts = set()
        for t in self.parsed:
            for e in t.sub_features:
                got_starts.add(int(e.location.start + 1))
        self.assertEquals(
            set([
                3038669,3186316,3190269,3192165,3192312,3192820,3195982,3195985,
                3196604,3197277,3199537,3203520,3203690,3204563,3215314,3335231,
                3338456,3343015,3344588,3345781,3347804,3349420,3359484,3367864,
                3369658,3382694,3385042,3411783,3456668,3503486,3638392,3648928,
                3660633,4280927,4333588,4341991,4342283,4350281,4399251,4481009,
                4483181,4483853,4485217,4486372,4763279,4767606,4772649,4774032,
                4775654
            ]),
            got_starts)

    def test_parsed_ends(self):
        got_ends = set()
        for t in self.parsed:
            for e in t.sub_features:
                got_ends.add(int(e.location.end))
        self.assertEquals(
            set([
                3038743,3186344,3190805,3192412,3193048,3197367,3197398,3199813,
                3205713,3206425,3207049,3215339,3335250,3338591,3343252,3344719,
                3345876,3347908,3349490,3359590,3368015,3369914,3382843,3385846,
                3411982,3456768,3503634,3640590,3648985,3661579,4283093,4340172,
                4342162,4342918,4350395,4399322,4482749,4483547,4483571,4483816,
                4483944,4486023,4486494,4764597,4766882,4767729,4772814,4774186,
                4775807,
            ]),
            got_ends)
        
