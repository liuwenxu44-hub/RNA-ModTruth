import itertools,math,unittest
from analyze_intervals import bin_geometry,ratio_extreme,envelopes,direction,decomposition_bounds

class IntervalTests(unittest.TestCase):
    def test_half_open_positive_endpoint(self):
        self.assertEqual(bin_geometry(191,.75)[4],[])
    def test_negative_threshold_singleton(self):
        lo,hi,mid,mandatory,parts=bin_geometry(64,.75)
        self.assertFalse(mid);self.assertFalse(mandatory);self.assertEqual(parts,[(.25,.25)])
        r=envelopes({(0,64):5},.75)
        self.assertEqual(r['minimum_S'],0);self.assertEqual(r['maximum_S'],5)
        self.assertEqual(r['variable_inf'],.25**2);self.assertFalse(r['variable_guaranteed_nonempty'])
    def test_baseline_retains_every_bin(self):
        self.assertTrue(all(bin_geometry(m,.5)[2] and bin_geometry(m,.5)[3] for m in range(256)))
    def test_ratio_mixed_optional_beats_all_in_or_out(self):
        self.assertAlmostEqual(ratio_extreme(1,.4,[(.1,1),(.9,1)]),.25)
        self.assertAlmostEqual(ratio_extreme(1,.4,[(.1,1),(.9,1)],True),.65)
    def test_ratio_optimum_matches_exhaustive_integer_counts(self):
        opts=[(.08,3),(.3,2),(.9,4)]
        values=[(.35+sum(k*x for k,(x,n) in zip(ks,opts)))/(1+sum(ks)) for ks in itertools.product(*[range(n+1) for x,n in opts])]
        self.assertAlmostEqual(ratio_extreme(1,.35,opts),min(values));self.assertAlmostEqual(ratio_extreme(1,.35,opts,True),max(values))
    def test_empty_and_potentially_empty(self):
        self.assertIsNone(envelopes({},.9)['variable_inf'])
        r=envelopes({(0,230):4},.9);self.assertFalse(r['variable_guaranteed_nonempty']);self.assertTrue(r['variable_possible_nonempty'])
        self.assertEqual(direction((r['variable_inf'],r['variable_sup']),(0.,1.),False)[0],'POTENTIALLY_UNDEFINED_REQUIRED_STRATUM')
    def test_fixed_membership_not_latent_selection(self):
        r=envelopes({(0,230):4},.9)
        self.assertLess(r['fixed_inf'],r['variable_inf']);self.assertAlmostEqual(r['variable_inf'],.9**2)
    def test_joint_membership_and_score_worlds_enclosed(self):
        h={(0,1):1,(0,230):2,(1,254):2,(1,25):1};r=envelopes(h,.9)
        for fs in itertools.product((0.,.5,.999),repeat=len(h)):
            points=[(y,(m+f)/256,n) for ((y,m),n),f in zip(h.items(),fs)]
            selected=[(y,q,n) for y,q,n in points if max(q,1-q)>=.9]
            loss=sum(n*(q-y)**2 for y,q,n in selected)/sum(n for y,q,n in selected)
            self.assertGreaterEqual(loss,r['variable_inf']-1e-12);self.assertLessEqual(loss,r['variable_sup']+1e-12)
    def test_overlap_is_not_reversal(self):
        self.assertEqual(direction((.1,.3),(.2,.4),True)[0],'NOT_CERTIFIED_BY_MARGINAL_ENVELOPES')
        self.assertEqual(direction((.1,.2),(.3,.4),True)[0],'CERTIFIED_INCREASE')
    def test_component_enclosures_and_baseline_identity(self):
        hs=[{(0,1):2,(0,230):2},{(1,254):3,(1,25):1}]
        base=[envelopes(h,.5) for h in hs];chosen=[envelopes(h,.9) for h in hs]
        for layer in ('fixed','variable'):
            r=decomposition_bounds(base,chosen,layer,.9)
            self.assertEqual(r['state'],'CONSERVATIVE_COMPONENT_ENCLOSURES')
            for f0,f1 in itertools.product((0.,.5,.999),repeat=2):
                points=[[(y,(m+f)/256,n) for (y,m),n in h.items()] for h,f in zip(hs,(f0,f1))]
                kept=[[(y,q,n) for y,q,n in pp if (max(q,1-q)>=.9 if layer=='variable' else max((int(q*256)+.5)/256,1-(int(q*256)+.5)/256)>=.9)] for pp in points]
                n0=[sum(n for y,q,n in pp) for pp in points];n1=[sum(n for y,q,n in pp) for pp in kept]
                l0=[sum(n*(q-y)**2 for y,q,n in pp)/n for pp,n in zip(points,n0)];l1=[sum(n*(q-y)**2 for y,q,n in pp)/n for pp,n in zip(kept,n1)]
                w0=[n/sum(n0) for n in n0];w1=[n/sum(n1) for n in n1]
                within=sum((a+b)/2*(d-c) for a,b,c,d in zip(w0,w1,l0,l1));weight=sum((c+d)/2*(b-a) for a,b,c,d in zip(w0,w1,l0,l1))
                self.assertLessEqual(r['within_inf']-1e-12,within);self.assertGreaterEqual(r['within_sup']+1e-12,within)
                self.assertLessEqual(r['weight_inf']-1e-12,weight);self.assertGreaterEqual(r['weight_sup']+1e-12,weight)
        self.assertEqual(decomposition_bounds(base,base,'variable',.5)['within_inf'],0.)

if __name__=='__main__':unittest.main(verbosity=2)
