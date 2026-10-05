import sys, unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
from robustness import Group

class DescriptiveMetrics(unittest.TestCase):
 def test_class_standardization(self):
  g=Group();g.add(0,.8,1);g.add(1,.9,9);r=g.row('fixture','fixture','canonical','fixture')
  self.assertAlmostEqual(r['Brier_observed'],.073)
  self.assertAlmostEqual(r['Brier_equal_class'],.325)
 def test_zero_positive_is_infinite(self):
  g=Group();g.add(0,0.,1);g.add(1,0.,1);r=g.row('x','x','missing_zero','x')
  self.assertEqual(r['log_loss'],'INFINITE');self.assertEqual(r['Brier_observed'],.5)
 def test_absent_class_not_zero(self):
  g=Group();g.add(1,.9,5);r=g.row('x','x','high_p_filter','x')
  self.assertIsNone(r['Brier_equal_class']);self.assertIsNone(r['FPR'])
 def test_interval_boundary(self):
  g=Group();g.add(0,.5/256,1,0.,1/256);g.add(1,.5/256,1,0.,1/256)
  r=g.row('x','x','canonical','x');self.assertEqual(r['log_loss_interval_supremum'],'INFINITE')
  self.assertLess(r['Brier_interval_infimum'],r['Brier_observed']);self.assertGreater(r['Brier_interval_supremum'],r['Brier_observed'])
 def test_site_mean_brier_identity(self):
  ps=[.1,.6,.8]; y=1;mu=sum(ps)/len(ps)
  raw=sum((p-y)**2 for p in ps)/len(ps); collapsed=(mu-y)**2
  self.assertAlmostEqual(raw-collapsed,sum((p-mu)**2 for p in ps)/len(ps))
 def test_ineligible_not_imputed_into_metrics(self):
  g=Group();g.den.update(opp=3,eligible=2,explicit=1);g.add(1,.8,1);g.add(1,0,1)
  r=g.row('x','x','missing_zero','x');self.assertEqual(r['scored'],2);self.assertEqual(r['opportunities'],3)
if __name__=='__main__':unittest.main(verbosity=2)
