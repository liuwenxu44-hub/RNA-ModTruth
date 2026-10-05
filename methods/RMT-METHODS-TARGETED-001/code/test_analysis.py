import unittest
from analyze import Group,retention,decompose,hiding,pooled

class AnalysisTests(unittest.TestCase):
    def test_half_open_upper_not_possible(self):
        self.assertEqual(retention(127,'one_sided_high_p_reference',.5),(False,False,False))
    def test_symmetric_baseline(self):
        self.assertTrue(all(retention(m,'symmetric_confidence',.5)==(True,True,True) for m in range(256)))
    def test_quantized_ambiguous(self):
        keep,d,p=retention(230,'symmetric_confidence',.9); self.assertTrue(keep); self.assertFalse(d); self.assertTrue(p)
    def test_complement_low_branch(self):
        self.assertEqual(retention(0,'symmetric_confidence',.99),(True,True,True))
    def test_decomposition(self):
        r=decompose({'a':(80,.2),'b':(20,.8)},{'a':(10,.3),'b':(90,.4)},{}); self.assertAlmostEqual(r[0]['addback_residual'],0)
    def test_disappearing_stratum(self):
        r=decompose({'a':(80,.2),'b':(20,.8)},{'a':(10,.3),'b':(0,None)},{})
        self.assertEqual(r[0]['state'],'UNIDENTIFIABLE'); self.assertEqual(r[1]['scope'],'COMMON_SUPPORT_ONLY'); self.assertEqual(r[1]['excluded_strata'],'b')
    def test_no_class_renormalization(self):
        a=Group(); a.opp=a.eligible=1; a.hist[0]=1; b=Group()
        self.assertIsNone(pooled([a.metric(0,'explicit_baseline',.5),b.metric(1,'explicit_baseline',.5)])['Brier_equal_class'])
    def test_controlled_hiding(self):
        g=Group(); g.hist.update({0:2,63:3,64:4,255:1})
        for y in (0,1):
            h=hiding(g,y); self.assertEqual(h['hidden'],5); self.assertLessEqual(h['lower_sum'],h['known_loss_sum']); self.assertGreaterEqual(h['upper_sum'],h['known_loss_sum'])
    def test_symmetric_t_not_monotonic_risk_assumption(self):
        g=Group(); g.hist.update({250:1,130:10})
        self.assertGreater(g.metric(0,'symmetric_confidence',.95)['brier'],g.metric(0,'explicit_baseline',.5)['brier'])

if __name__=='__main__': unittest.main(verbosity=2)
