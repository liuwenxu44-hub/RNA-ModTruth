"""Synthetic syntax/semantics positives and negatives; no biological gold claims."""
import array
import copy
import json
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'code'))
import pysam
from canary import decode_mm,cigar_pairs
from metrics import calculate


def read(seq,mm,ml=None,flag=0):
    item=pysam.AlignedSegment()
    item.query_name='synthetic-parser-fixture'
    item.query_sequence=seq
    item.flag=flag
    item.reference_id=0
    item.reference_start=10
    item.cigartuples=[(0,len(seq))]
    item.set_tag('MM',mm)
    if ml is not None:
        item.set_tag('ML',array.array('B',ml))
    return item


class ParserTests(unittest.TestCase):
    def test_explicit_delta(self):
        d,_,_=decode_mm(read('AACAA','A+a.,0,1;',[0,255]))
        self.assertEqual(d,{(0,'A',0,'a'):0,(3,'A',0,'a'):255})

    def test_empty_numeric_group_consumes_no_ml(self):
        d,m,_=decode_mm(read('ATACA','A+a.,0;T+17802.;',[17]))
        self.assertEqual(d,{(0,'A',0,'a'):17})
        self.assertEqual(m[('T',0,17802)],'.')

    def test_question_dot_omitted_distinct(self):
        for header,expected in [('A+a?;','?'),('A+a.;','.'),('A+a;','omitted')]:
            d,m,_=decode_mm(read('AAA',header,[]))
            self.assertEqual(d,{})
            self.assertEqual(m[('A',0,'a')],expected)

    def test_missing_ml_is_not_zero(self):
        d,_,_=decode_mm(read('AAA','A+a.,0;'))
        self.assertIsNone(d[(0,'A',0,'a')])

    def test_reverse_original_orientation(self):
        item=read('TTGTT','A+a.,1;',[129],16)
        d,_,_=decode_mm(item)
        self.assertEqual(d,{(3,'A',1,'a'):129})
        native={(q,b,s,c):v for (b,s,c),pairs in item.modified_bases.items() for q,v in pairs}
        self.assertEqual(d,native)

    def test_interleaved_multicode(self):
        d,_,_=decode_mm(read('AAA','A+ah.,0,1;',[10,20,30,40]))
        self.assertEqual(d,{(0,'A',0,'a'):10,(0,'A',0,'h'):20,(2,'A',0,'a'):30,(2,'A',0,'h'):40})

    def test_stale_mn(self):
        item=read('AAA','A+a.,0;',[128])
        item.set_tag('MN',4)
        with self.assertRaisesRegex(ValueError,'STALE_MN'):
            decode_mm(item)

    def test_bad_delta(self):
        with self.assertRaisesRegex(ValueError,'MM_OFFSET'):
            decode_mm(read('AAA','A+a.,3;',[128]))

    def test_ml_length(self):
        with self.assertRaisesRegex(ValueError,'ML_TOO_SHORT'):
            decode_mm(read('AAA','A+a.,0,0;',[128]))
        with self.assertRaisesRegex(ValueError,'ML_LENGTH'):
            decode_mm(read('AAA','A+a.,0;',[128,129]))

    def test_cigar_gaps_and_clips(self):
        item=read('AACCGG','A+a.,0;',[128])
        item.cigartuples=[(4,1),(0,2),(2,1),(1,1),(0,2)]
        self.assertEqual(cigar_pairs(item),item.get_aligned_pairs())

    def test_no_tag(self):
        item=read('AAA','A+a.,0;',[128])
        item.set_tag('MM',None)
        d,m,state=decode_mm(item)
        self.assertEqual((d,m,state),({},{},'NO_MM_TAG'))


if __name__=='__main__':
    unittest.main(verbosity=2)
