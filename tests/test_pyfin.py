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

    def test_binomial_tree_discrete_dividend_valuation(self):
        mat = 9.0 / 12
        num_steps = 9 * 4
        self.assertAlmostEquals(13.93, Option(opt_type=OptionType.PUT, spot0=90, strike=100, mat=mat,
                                              vol=0.35, riskless_rate=0.10) \
                                .run_model(model=OptionModel.BINOMIAL_TREE, num_steps=num_steps)[OptionMeasure.VALUE], delta=0.01)
        self.assertAlmostEquals(13.93, Option(opt_type=OptionType.PUT, spot0=90, strike=100, mat=mat,
                                              vol=0.35, riskless_rate=0.10, yield_=lambda t1, t2: 0) \
                                .run_model(model=OptionModel.BINOMIAL_TREE, num_steps=num_steps)[OptionMeasure.VALUE], delta=0.01)
        self.assertAlmostEquals(14.77, Option(opt_type=OptionType.PUT, spot0=90, strike=100, mat=mat,
                                              vol=0.35, riskless_rate=0.10, yield_=0.04) \
                                .run_model(model=OptionModel.BINOMIAL_TREE, num_steps=num_steps)[OptionMeasure.VALUE], delta=0.01)

        def div1_fn(t1, t2):  # $2.50 dividend after about a one month
            n = t2 / (mat / num_steps)
            return 2.50 if n == 4 else 0
        self.assertAlmostEquals(15.45, Option(opt_type=OptionType.PUT, spot0=90, strike=100, mat=mat,
                                              vol=0.35, riskless_rate=0.10, yield_=div1_fn) \
                                .run_model(model=OptionModel.BINOMIAL_TREE, num_steps=num_steps)[OptionMeasure.VALUE], delta=0.01)

        def div_monthly_fn(t1, t2):  # $2.50 dividend every month
            n = t2 / (mat / num_steps)
            return 2.50 if (n % 4) == 0 else 0
        self.assertAlmostEquals(3.57, Option(opt_type=OptionType.CALL, spot0=90, strike=100, mat=mat,
                                             vol=0.35, riskless_rate=0.10, yield_=div_monthly_fn) \
                                .run_model(model=OptionModel.BINOMIAL_TREE, num_steps=num_steps)[OptionMeasure.VALUE], delta=0.01)
        self.assertAlmostEquals(24.51, Option(opt_type=OptionType.PUT, spot0=90, strike=100, mat=mat,
                                              vol=0.35, riskless_rate=0.10, yield_=div_monthly_fn) \
                                .run_model(model=OptionModel.BINOMIAL_TREE, num_steps=num_steps)[OptionMeasure.VALUE], delta=0.01)

    def test_monte_carlo_valuation(self):
        self.assertAlmostEquals(2.51, Option(opt_type=OptionType.CALL, spot0=90, strike=100, mat=3.0 / 12,
                                             vol=0.35, riskless_rate=0.10, yield_=0.05, exer_type=OptionExerciseType.EUROPEAN) \
                                .run_model(model=OptionModel.MONTE_CARLO, random_seed=12345)[OptionMeasure.VALUE], delta=0.01)
        self.assertAlmostEquals(3.60, Option(opt_type=OptionType.CALL, spot0=90, strike=100, mat=3.0 / 12,
                                             vol=0.35, riskless_rate=0.10, yield_=0.05, exer_type=OptionExerciseType.AMERICAN) \
                                .run_model(model=OptionModel.MONTE_CARLO, random_seed=12345)[OptionMeasure.VALUE], delta=0.01)

        mat = 9.0 / 12
        num_steps = 9 * 4
        def div_monthly_fn(t1, t2):  # $2.50 dividend every month
            n = t2 / (mat / num_steps)
            return 2.50 if (n % 4) == 0 else 0

        self.assertAlmostEquals(10.16, Option(opt_type=OptionType.PUT, spot0=90, strike=100, mat=mat,
                                              vol=0.35, riskless_rate=0.10, exer_type=OptionExerciseType.EUROPEAN) \
                                .run_model(model=OptionModel.MONTE_CARLO, num_steps=num_steps, random_seed=12345)[OptionMeasure.VALUE], delta=0.01)
        self.assertAlmostEquals(23.39, Option(opt_type=OptionType.PUT, spot0=90, strike=100, mat=mat,
                                              vol=0.35, riskless_rate=0.10, yield_=div_monthly_fn, exer_type=OptionExerciseType.EUROPEAN) \
                                .run_model(model=OptionModel.MONTE_CARLO, num_steps=num_steps, random_seed=12345)[OptionMeasure.VALUE], delta=0.01)
        self.assertAlmostEquals(1.40, Option(opt_type=OptionType.CALL, spot0=90, strike=100, mat=mat,
                                             vol=0.35, riskless_rate=0.10, yield_=div_monthly_fn, exer_type=OptionExerciseType.EUROPEAN) \
                                .run_model(model=OptionModel.MONTE_CARLO, num_steps=num_steps, random_seed=12345)[OptionMeasure.VALUE], delta=0.01)

        greek_put_opt_results = Option(opt_type=OptionType.PUT, spot0=90, strike=100, mat=3.0 / 12,
                                       vol=0.35, riskless_rate=0.10, yield_=0.05, exer_type=OptionExerciseType.EUROPEAN) \
            .run_model(model=OptionModel.MONTE_CARLO, random_seed=12345)
        self.assertAlmostEqual(-0.6442, greek_put_opt_results[OptionMeasure.DELTA], delta=0.0001)
        self.assertAlmostEqual(-6.7226, greek_put_opt_results[OptionMeasure.THETA], delta=0.0001)
        self.assertAlmostEqual(-17.0219, greek_put_opt_results[OptionMeasure.RHO], delta=0.0001)
        self.assertAlmostEqual(+15.2718, greek_put_opt_results[OptionMeasure.VEGA], delta=0.0001)
        self.assertAlmostEqual(+0.0254, greek_put_opt_results[OptionMeasure.GAMMA], delta=0.0001)

        greek_call_opt_results = Option(opt_type=OptionType.CALL, spot0=90, strike=100, mat=3.0 / 12,
                                       vol=0.35, riskless_rate=0.10, yield_=0.05, exer_type=OptionExerciseType.AMERICAN) \
            .run_model(model=OptionModel.MONTE_CARLO, random_seed=12345)
        self.assertAlmostEqual(+0.4217, greek_call_opt_results[OptionMeasure.DELTA], delta=0.0001)
        self.assertAlmostEqual(+60.4526, greek_call_opt_results[OptionMeasure.THETA], delta=0.0001)
        self.assertAlmostEqual(+6.4028, greek_call_opt_results[OptionMeasure.RHO], delta=0.0001)
        self.assertAlmostEqual(-55.4468, greek_call_opt_results[OptionMeasure.VEGA], delta=0.0001)
        self.assertAlmostEqual(+0.0163, greek_call_opt_results[OptionMeasure.GAMMA], delta=0.0001)

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

    def test_model_cache(self):
        opt = Option(opt_type=OptionType.CALL, spot0=90, strike=100, mat=3.0 / 12,
                     vol=0.35, riskless_rate=0.10, yield_=0.025)

        self.assertEquals(0, len(opt.model_cache))
        self.assertAlmostEquals(3.19, opt.run_model(model=OptionModel.BINOMIAL_TREE)[OptionMeasure.VALUE], delta=0.01)
        self.assertEquals([OptionModel.BINOMIAL_TREE], opt.model_cache.keys())
        self.assertAlmostEquals(3.19, opt.run_model(model=OptionModel.BINOMIAL_TREE)[OptionMeasure.VALUE], delta=0.01)

        opt.model_cache[OptionModel.BINOMIAL_TREE][OptionMeasure.VALUE] = 42.42
        self.assertAlmostEquals(42.42, opt.run_model(model=OptionModel.BINOMIAL_TREE)[OptionMeasure.VALUE], delta=0.01)

        opt.vol = 0.55
        self.assertAlmostEquals(6.76, opt.run_model(model=OptionModel.BINOMIAL_TREE)[OptionMeasure.VALUE], delta=0.01)

    def tearDown(self):
        pass

if __name__ == '__main__':
    unittest.main()
