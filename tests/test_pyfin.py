#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
test_pyfin
----------------------------------

Unit tests the valuation & risk sensitivities for options in the `pyfin` module.
"""

import unittest

from pyfin.pyfin import OptionType, OptionExerciseType, Option, OptionModel, OptionMeasure


class TestPyfin(unittest.TestCase):

    def setUp(self):
        pass

    def test_black_scholes_valuation(self):
        self.assertAlmostEquals(13.70, Option(opt_type=OptionType.CALL, spot0=100, strike=95, mat=3.0 / 12,
                                              vol=0.50, riskless_rate=0.10, exer_type=OptionExerciseType.EUROPEAN).run_model(model=OptionModel.BLACK_SCHOLES)[OptionMeasure.VALUE], delta=0.01)
        self.assertAlmostEquals(10.03, Option(opt_type=OptionType.PUT, spot0=90, strike=100, mat=1.0 / 12,
                                              vol=0.35, riskless_rate=0.10, exer_type=OptionExerciseType.EUROPEAN).run_model(model=OptionModel.BLACK_SCHOLES)[OptionMeasure.VALUE], delta=0.01)
        self.assertAlmostEquals(3.45, Option(opt_type=OptionType.CALL, spot0=90, strike=100, mat=3.0 / 12,
                                             vol=0.35, riskless_rate=0.10, exer_type=OptionExerciseType.EUROPEAN).run_model(model=OptionModel.BLACK_SCHOLES)[OptionMeasure.VALUE], delta=0.01)
        self.assertAlmostEquals(3.06, Option(opt_type=OptionType.CALL, spot0=90, strike=100, mat=3.0 / 12,
                                             vol=0.35, riskless_rate=0.10, yield_=0.05, exer_type=OptionExerciseType.EUROPEAN).run_model(model=OptionModel.BLACK_SCHOLES)[OptionMeasure.VALUE], delta=0.01)

        greek_call_opt = Option(opt_type=OptionType.CALL, spot0=90, strike=100, mat=3.0 / 12,
                                vol=0.35, riskless_rate=0.10, yield_=0.05, exer_type=OptionExerciseType.EUROPEAN)
        greek_put_opt = Option(opt_type=OptionType.PUT, spot0=90, strike=100, mat=3.0 / 12,
                               vol=0.35, riskless_rate=0.10, yield_=0.05, exer_type=OptionExerciseType.EUROPEAN)
        self.assertAlmostEquals(+0.3248,   greek_call_opt.run_model(model=OptionModel.BLACK_SCHOLES)[OptionMeasure.DELTA], delta=0.0001)
        self.assertAlmostEquals(-0.6628,   greek_put_opt.run_model(model=OptionModel.BLACK_SCHOLES) [OptionMeasure.DELTA], delta=0.0001)
        self.assertAlmostEquals(-12.4048,  greek_call_opt.run_model(model=OptionModel.BLACK_SCHOLES)[OptionMeasure.THETA], delta=0.0001)
        self.assertAlmostEquals(-7.0958,   greek_put_opt.run_model(model=OptionModel.BLACK_SCHOLES) [OptionMeasure.THETA], delta=0.0001)
        self.assertAlmostEquals(+6.5405,   greek_call_opt.run_model(model=OptionModel.BLACK_SCHOLES)[OptionMeasure.RHO],   delta=0.0001)
        self.assertAlmostEquals(-17.8422,  greek_put_opt.run_model(model=OptionModel.BLACK_SCHOLES) [OptionMeasure.RHO],   delta=0.0001)
        self.assertAlmostEquals(+16.0714,  greek_call_opt.run_model(model=OptionModel.BLACK_SCHOLES)[OptionMeasure.VEGA],  delta=0.0001)
        self.assertAlmostEquals(+16.0714,  greek_put_opt.run_model(model=OptionModel.BLACK_SCHOLES) [OptionMeasure.VEGA],  delta=0.0001)
        self.assertAlmostEquals(+0.0227,   greek_call_opt.run_model(model=OptionModel.BLACK_SCHOLES)[OptionMeasure.GAMMA], delta=0.0001)
        self.assertAlmostEquals(+0.0227,   greek_put_opt.run_model(model=OptionModel.BLACK_SCHOLES) [OptionMeasure.GAMMA], delta=0.0001)

    def test_binomial_tree_valuation(self):
        self.assertAlmostEquals(10.32, Option(opt_type=OptionType.PUT, spot0=90, strike=100, mat=1.0 / 12,
                                              vol=0.35, riskless_rate=0.10).run_model(model=OptionModel.BINOMIAL_TREE)[OptionMeasure.VALUE], delta=0.01)
        self.assertAlmostEquals(10.01, Option(opt_type=OptionType.PUT, spot0=90, strike=100, mat=1.0 / 12,
                                              vol=0.35, riskless_rate=0.10, exer_type=OptionExerciseType.EUROPEAN).run_model(model=OptionModel.BINOMIAL_TREE)[OptionMeasure.VALUE], delta=0.01)
        self.assertAlmostEquals(10.40, Option(opt_type=OptionType.PUT, spot0=90, strike=100, mat=1.0 / 12,
                                              vol=0.35, riskless_rate=0.10, yield_=0.025).run_model(model=OptionModel.BINOMIAL_TREE)[OptionMeasure.VALUE], delta=0.01)
        self.assertAlmostEquals(3.19, Option(opt_type=OptionType.CALL, spot0=90, strike=100, mat=3.0 / 12,
                                             vol=0.35, riskless_rate=0.10, yield_=0.025).run_model(model=OptionModel.BINOMIAL_TREE)[OptionMeasure.VALUE], delta=0.01)

        greek_call_opt = Option(opt_type=OptionType.CALL, spot0=90, strike=100, mat=3.0 / 12,
                                vol=0.35, riskless_rate=0.10, yield_=0.05, exer_type=OptionExerciseType.AMERICAN)
        greek_put_opt = Option(opt_type=OptionType.PUT, spot0=90, strike=100, mat=3.0 / 12,
                               vol=0.35, riskless_rate=0.10, yield_=0.05, exer_type=OptionExerciseType.AMERICAN)
        self.assertAlmostEquals(+0.3205,   greek_call_opt.run_model(model=OptionModel.BINOMIAL_TREE)[OptionMeasure.DELTA], delta=0.0001)
        self.assertAlmostEquals(-0.6973,   greek_put_opt.run_model(model=OptionModel.BINOMIAL_TREE) [OptionMeasure.DELTA], delta=0.0001)
        self.assertAlmostEquals(-12.6706,  greek_call_opt.run_model(model=OptionModel.BINOMIAL_TREE)[OptionMeasure.THETA], delta=0.0001)
        self.assertAlmostEquals(-8.3641,   greek_put_opt.run_model(model=OptionModel.BINOMIAL_TREE) [OptionMeasure.THETA], delta=0.0001)
        self.assertAlmostEquals(+6.4561,   greek_call_opt.run_model(model=OptionModel.BINOMIAL_TREE)[OptionMeasure.RHO],   delta=0.0001)
        self.assertAlmostEquals(-10.9031,  greek_put_opt.run_model(model=OptionModel.BINOMIAL_TREE) [OptionMeasure.RHO],   delta=0.0001)
        self.assertAlmostEquals(+14.1128,  greek_call_opt.run_model(model=OptionModel.BINOMIAL_TREE)[OptionMeasure.VEGA],  delta=0.0001)
        self.assertAlmostEquals(+13.7077,  greek_put_opt.run_model(model=OptionModel.BINOMIAL_TREE) [OptionMeasure.VEGA],  delta=0.0001)
        self.assertAlmostEquals(+0.0268,   greek_call_opt.run_model(model=OptionModel.BINOMIAL_TREE)[OptionMeasure.GAMMA], delta=0.0001)
        self.assertAlmostEquals(+0.0178,   greek_put_opt.run_model(model=OptionModel.BINOMIAL_TREE) [OptionMeasure.GAMMA], delta=0.0001)


    def test_implied_volatility(self):
        self.assertAlmostEquals(0.30, Option.imply_volatility(10.10, opt_type=OptionType.PUT, spot0=90, strike=100, mat=1.0 / 12,
                                                              riskless_rate=0.10), delta=0.01)
        self.assertAlmostEquals(0.52, Option.imply_volatility(2.25, opt_type=OptionType.CALL, spot0=90, strike=100, mat=1.0 / 12,
                                                              riskless_rate=0.10), delta=0.01)

        premium = 10.10
        opt_kwargs = {
            'opt_type': OptionType.PUT,
            'spot0': 90,
            'strike': 100,
            'mat': 1.0 / 12,
            'riskless_rate': 0.10,
        }
        self.assertAlmostEquals(0.30, Option.imply_volatility(premium, **opt_kwargs), delta=0.01)
        opt_kwargs['vol'] = Option.imply_volatility(premium, **opt_kwargs)
        self.assertAlmostEquals(premium, Option(**opt_kwargs).run_model(model=OptionModel.BINOMIAL_TREE)[OptionMeasure.VALUE], delta=0.01)


    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
