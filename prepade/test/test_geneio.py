import unittest
import numpy as np
from prepade.geneio import read_exons, ExonIndex

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
            [set([('chr1',  4, 10)]),
             set([('chr1', 10, 20)]),
             set([('chr1', 10, 20), ('chr1', 14, 30)]),
             set([('chr1', 10, 20), ('chr1', 14, 30), ('chr1', 17, 25)]),
             set([('chr1', 14, 30), ('chr1', 17, 25)]),
             set([('chr1', 14, 30)]),
             set(),
             set([('chr1', 34, 40)]),
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
