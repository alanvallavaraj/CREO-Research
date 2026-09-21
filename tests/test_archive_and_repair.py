"""Unit tests of the ACTUAL archived repair functions and archived data integrity.

Uses AST to extract the three unmodified repair functions without executing the
experiment module's import-time external-CEC extraction or 30x29x10000D run loop.
"""
import ast
import unittest
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'ablation_experiment.py'
NAMES={'clip_to_bounds','total_violation','repair_solution'}

def load_actual_repair():
    tree=ast.parse(SOURCE.read_text(encoding='utf-8'))
    nodes=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in NAMES]
    if {n.name for n in nodes}!=NAMES:raise ValueError('Repair functions missing or renamed')
    env={'np':np,'LOWER_BOUND':-100.0,'UPPER_BOUND':100.0}
    exec(compile(ast.Module(body=nodes,type_ignores=[]),str(SOURCE),'exec'),env)
    return env

class OriginalRepairTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.m=load_actual_repair()
    def test_clip_bounds(self):
        self.assertTrue(np.array_equal(self.m['clip_to_bounds'](np.array([-120.,12.,150.])),np.array([-100.,12.,100.])))
    def test_identity_preferred_region(self):
        x=np.array([-30.,22.,51.]);self.assertEqual(self.m['total_violation'](x),0.)
        self.assertTrue(np.array_equal(self.m['repair_solution'](x),x))
    def test_shrink_and_displacement(self):
        x=np.array([90.,90.]);y=self.m['repair_solution'](x)
        self.assertTrue(np.allclose(y,np.array([57.6,57.6])))
        self.assertTrue(np.all(y-x<0))
        self.assertEqual(self.m['total_violation'](y),0.)
    def test_repair_is_feasible_after_three_steps(self):
        for d in (10,30):
            for x in [np.full(d,100.),np.full(d,-100.),np.linspace(-100,100,d)]:
                self.assertLessEqual(self.m['total_violation'](self.m['repair_solution'](x)),1e-12)
    def test_archive_shape(self):
        import pandas as pd
        for d in (10,30):
            m=pd.read_csv(ROOT/f'cec17_results_D{d}/all_benchmark_results.csv')
            a=pd.read_csv(ROOT/f'cec17_ablation_small_D{d}/all_ablation_results.csv')
            self.assertEqual(len(m),29*7);self.assertEqual(len(a),29*3)
            p=a.pivot(index='Function',columns='Algorithm',values='Mean')
            w=(p['CREO-DE']<p['CREO-DE-NoRG']).sum()
            self.assertEqual(int(w),12 if d==10 else 1)

if __name__=='__main__':unittest.main()
