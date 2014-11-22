# -*- coding: utf-8 -*-
import abc
from math import exp
import math
import numbers
import time

import numpy as np
import scipy.optimize
import scipy.stats


def enum(**enums):
    return type('Enum', (), enums)

OptionType = enum(CALL='call', PUT='put')
OptionExerciseType = enum(EUROPEAN='european', AMERICAN='american')
OptionModel = enum(BLACK_SCHOLES='black_scholes', BINOMIAL_TREE='binomial_tree', MONTE_CARLO='monte_carlo')
OptionMeasure = enum(VALUE='value', DELTA='delta', THETA='theta', RHO='rho', VEGA='vega', GAMMA='gamma')

DEFAULT_BINOMIAL_TREE_NUM_STEPS = 25
DEFAULT_MONTE_CARLO_NUM_STEPS = 50
DEFAULT_MONTE_CARLO_NUM_PATHS = 100


class Instrument(object):
    @abc.abstractmethod
    def run_model(self):
        """Calculate measures (i.e. theoretical value & greeks) for this instrument"""


class Option(Instrument):
    def __init__(self, opt_type, spot0, strike, mat, vol=None, riskless_rate=None, yield_=None, exer_type=None):
        self.opt_type = opt_type
        self.spot0 = spot0
        self.strike = strike
        self.vol = vol
        self.mat = mat
        self.riskless_rate = riskless_rate or 0
        self.yield_ = yield_ or 0
        self.exer_type = exer_type or OptionExerciseType.AMERICAN

        self.model_cache = {}
        self.model_cache_param_hashes = {}

    def copy(self):
        return Option(self.opt_type, self.spot0, self.strike, self.mat,
                      vol=self.vol, riskless_rate=self.riskless_rate, yield_=self.yield_, exer_type=self.exer_type)

    def param_hash(self):
        return hash((self.opt_type, self.spot0, self.strike, self.mat,
                     self.vol, self.riskless_rate, self.yield_, self.exer_type))

    @staticmethod
    def imply_volatility(premium, max_vol=0.99, *args, **kwargs):
        def obj_fn(vol_guess):
            kwargs['vol'] = vol_guess
            val = Option(*args, **kwargs).run_model()[OptionMeasure.VALUE]
            return val - premium
        try:
            return scipy.optimize.bisect(obj_fn, 0.01, max_vol, xtol=0.0025)
        except ValueError:
            return None

    def run_model(self, model=OptionModel.BINOMIAL_TREE, **kwargs):
        curr_param_hash = self.param_hash()
        if model in self.model_cache:
            prev_param_hash = self.model_cache_param_hashes[model]
            if self.param_hash() == prev_param_hash:
                return self.model_cache[model]

        if model == OptionModel.BLACK_SCHOLES:
            result = self.black_scholes(**kwargs)
        elif model == OptionModel.BINOMIAL_TREE:
            result = self.binomial_tree(**kwargs)
        else:  # model == OptionModel.MONTE_CARLO:
            result = self.monte_carlo(**kwargs)

        self.model_cache[model] = result
        self.model_cache_param_hashes[model] = curr_param_hash
        return result

    def black_scholes(self):
        assert self.exer_type == OptionExerciseType.EUROPEAN, \
            "Black-Scholes does not support early exercise"
        assert (self.yield_ is None) or is_number(self.yield_), \
            "Black-Scholes does not support discrete dividends"

        sqrt_mat = self.mat ** 0.5
        d1 = (math.log(float(self.spot0) / self.strike) + (self.riskless_rate - self.yield_ + 0.5 * self.vol * self.vol) * self.mat) / (self.vol * sqrt_mat)
        d2 = d1 - self.vol * (self.mat ** 0.5)
        d1_pdf = scipy.stats.norm.pdf(d1)
        riskless_disc = exp(-self.riskless_rate * self.mat)
        yield_disc = exp(-self.yield_ * self.mat)
        if self.opt_type == OptionType.CALL:
            d1_cdf = scipy.stats.norm.cdf(d1)
            d2_cdf = scipy.stats.norm.cdf(d2)
            delta = yield_disc * d1_cdf
            val = self.spot0 * delta - riskless_disc * self.strike * d2_cdf
            theta = -yield_disc * (self.spot0 * d1_pdf * self.vol) / (2 * sqrt_mat) - \
                self.riskless_rate * self.strike * riskless_disc * d2_cdf + \
                self.yield_ * self.spot0 * yield_disc * d1_cdf
            rho = self.strike * self.mat * riskless_disc * d2_cdf

        else:  # opt_type == OptionType.PUT:
            neg_d1_cdf = scipy.stats.norm.cdf(-d1)
            neg_d2_cdf = scipy.stats.norm.cdf(-d2)
            delta = -yield_disc * neg_d1_cdf
            val = riskless_disc * self.strike * neg_d2_cdf + self.spot0 * delta
            theta = -yield_disc * (self.spot0 * d1_pdf * self.vol) / (2 * sqrt_mat) + \
                self.riskless_rate * self.strike * riskless_disc * neg_d2_cdf - \
                self.yield_ * self.spot0 * yield_disc * neg_d1_cdf
            rho = -self.strike * self.mat * riskless_disc * neg_d2_cdf

        vega = self.spot0 * yield_disc * d1_pdf * sqrt_mat
        gamma = yield_disc * (d1_pdf / (self.spot0 * self.vol * sqrt_mat))

        return {
            OptionMeasure.VALUE: val,
            OptionMeasure.DELTA: delta,
            OptionMeasure.THETA: theta,
            OptionMeasure.RHO: rho,
            OptionMeasure.VEGA: vega,
            OptionMeasure.GAMMA: gamma,
        }

    def binomial_tree(self, num_steps=None, sens_degree=2):
        val_num_steps = num_steps or DEFAULT_BINOMIAL_TREE_NUM_STEPS
        val_cache = {}

        dt = float(self.mat) / val_num_steps
        up_fact = exp(self.vol * (dt ** 0.5))
        down_fact = 1.0 / up_fact
        cont_yield = 0 if hasattr(self.yield_, '__call__') else self.yield_
        prob_up = (exp((self.riskless_rate - cont_yield) * dt) - down_fact) / (up_fact - down_fact)

        def divs_pv(n):
            if not hasattr(self.yield_, '__call__'): return 0
            pv_divs = 0
            for step_num in xrange(n, val_num_steps):
                div_t1, div_t2 = dt * step_num, dt * (step_num + 1)
                div = self.yield_(div_t1, div_t2)
                if is_number(div):
                    mid_div_t = 0.5 * (div_t1 + div_t2)
                    pv_divs += exp(-self.riskless_rate * mid_div_t) * div
            return pv_divs

        bin_spot0 = self.spot0 - divs_pv(0)

        def spot_price(n, num_ups, num_downs):
            return (bin_spot0 + divs_pv(n)) * (up_fact ** (num_ups - num_downs))

        def node_value(n, num_ups, num_downs):
            value_cache_key = (n, num_ups, num_downs)
            if value_cache_key not in val_cache:
                spot = spot_price(n, num_ups, num_downs)
                if self.opt_type == OptionType.CALL:
                    exer_profit = max(0, spot - self.strike)
                else:  # opt_type == OptionType.PUT:
                    exer_profit = max(0, self.strike - spot)

                if n >= val_num_steps:
                    val = exer_profit
                else:
                    fv = prob_up * node_value(n + 1, num_ups + 1, num_downs) + \
                        (1 - prob_up) * node_value(n + 1, num_ups, num_downs + 1)
                    pv = exp(-self.riskless_rate * dt) * fv
                    if self.exer_type == OptionExerciseType.AMERICAN:
                        val = max(pv, exer_profit)
                    else:  # exer_type == OptionExerciseType.EUROPEAN
                        val = pv

                val_cache[value_cache_key] = val
            return val_cache[value_cache_key]

        val = node_value(0, 0, 0)
        delta, theta, rho, vega, gamma = None, None, None, None, None
        if sens_degree >= 1:
            delta = (node_value(1, 1, 0) - node_value(1, 0, 1)) / (bin_spot0 * up_fact - bin_spot0 * down_fact)
            theta = (node_value(2, 1, 1) - val) / (2 * dt)

            rho = self.sensitivity('riskless_rate', 0.0001,
                                   model=OptionModel.BINOMIAL_TREE, sens_degree=sens_degree - 1)
            vega = self.sensitivity('vol', 0.001,
                                    model=OptionModel.BINOMIAL_TREE, sens_degree=sens_degree - 1)

            delta_up = (node_value(2, 2, 0) - node_value(2, 1, 1)) / (bin_spot0 * up_fact - bin_spot0 * down_fact)
            delta_down = (node_value(2, 1, 1) - node_value(2, 0, 2)) / (bin_spot0 * up_fact - bin_spot0 * down_fact)
            gamma = (delta_up - delta_down) / ((bin_spot0 * up_fact * up_fact - bin_spot0 * down_fact * down_fact) / 2)

        return {
            OptionMeasure.VALUE: val,
            OptionMeasure.DELTA: delta,
            OptionMeasure.THETA: theta,
            OptionMeasure.RHO: rho,
            OptionMeasure.VEGA: vega,
            OptionMeasure.GAMMA: gamma,
        }

    def monte_carlo(self, num_steps=None, num_paths=None, random_seed=None, sens_degree=2):
        val_num_steps = num_steps or DEFAULT_MONTE_CARLO_NUM_STEPS
        val_num_paths = num_paths or DEFAULT_MONTE_CARLO_NUM_PATHS
        val_random_seed = random_seed or int(time.time() * 1000)

        np.random.seed(val_random_seed)
        dt = self.mat / val_num_steps
        sqrt_dt = dt ** 0.5
        cont_yield = 0 if hasattr(self.yield_, '__call__') else self.yield_
        drift = (self.riskless_rate - cont_yield - (self.vol ** 2) / 2) * dt

        def cast_spot_paths():
            result = np.empty([val_num_paths, val_num_steps])
            for path_num in xrange(val_num_paths):
                result1 = np.empty([val_num_steps])
                spot = self.spot0
                result1[0] = spot
                for step_num in xrange(1, val_num_steps):
                    spot *= exp(drift + self.vol * np.random.normal() * sqrt_dt)
                    if hasattr(self.yield_, '__call__'):
                        div_t1, div_t2 = dt * step_num, dt * (step_num + 1)
                        div = self.yield_(div_t1, div_t2)
                        if is_number(div): spot -= div
                    result1[step_num] = spot

                result[path_num] = result1
            return result

        def step_exercise_profits(spot_paths, step_num):
            exer_profits = np.empty([val_num_paths])
            for path_num in xrange(val_num_paths):
                spot = spot_paths[path_num, step_num]
                if self.opt_type == OptionType.CALL:
                    exer_profit = max(0, spot - self.strike)
                else:  # opt_type == OptionType.PUT:
                    exer_profit = max(0, self.strike - spot)
                exer_profits[path_num] = exer_profit
            return exer_profits

        if self.exer_type == OptionExerciseType.EUROPEAN:
            spot_paths = cast_spot_paths()
            exer_profits = step_exercise_profits(spot_paths, -1)
            pvs = exp(-self.riskless_rate * self.mat) * exer_profits
            val = pvs.mean()

        else:
            spot_paths = cast_spot_paths()
            cashflow_paths = np.zeros([val_num_paths, val_num_steps])
            cashflow_paths[:, val_num_steps - 1] = step_exercise_profits(spot_paths, -1)

            def cashflow_path_pv(path_num, from_step_num=0):
                from_cashflow_path = cashflow_paths[path_num, from_step_num:]
                max_step_num = from_cashflow_path.argmax()
                max_cashflow = from_cashflow_path[max_step_num]
                return exp(-self.riskless_rate * (dt * max_step_num)) * max_cashflow

            step_models = [None] * val_num_steps
            for step_num in xrange(val_num_steps - 2, -1, -1):
                exer_profits = step_exercise_profits(spot_paths, step_num)
                xs, ys = [], []
                for path_num in xrange(val_num_paths):
                    if exer_profits[path_num] <= 0: continue
                    xs += [spot_paths[path_num, step_num]]
                    cont_val = cashflow_path_pv(path_num, step_num + 1)
                    ys += [exp(-self.riskless_rate * dt) * cont_val]
                if len(xs) == 0: continue

                basis_xs = [xs, np.multiply(xs, xs), [1] * len(xs)]
                [x_coeff, x_sq_coeff, int_coeff] = np.linalg.lstsq(np.transpose(np.array(basis_xs)), ys)[0]
                step_models[step_num] = [x_coeff, x_sq_coeff, int_coeff]

                for path_num in xrange(val_num_paths):
                    exer_profit = exer_profits[path_num]
                    if exer_profit <= 0: continue
                    spot = spot_paths[path_num, step_num]

                    cont_pv = (x_coeff * spot) + (x_sq_coeff * spot * spot) + int_coeff
                    if exer_profit > cont_pv:
                        cashflow_paths[path_num, step_num] = exer_profit
                        cashflow_paths[path_num, (step_num + 1):] = 0

            pvs = []
            for path_num in xrange(val_num_paths):
                step_num = cashflow_paths[path_num].argmax()
                cashflow = cashflow_paths[path_num][step_num]
                pvs += [exp(-self.riskless_rate * (dt * step_num)) * cashflow]
            val = np.array(pvs).mean()

        delta, theta, rho, vega, gamma = None, None, None, None, None
        if sens_degree >= 1:
            delta = self.sensitivity('spot0', self.spot0 * 1.0e-05, model=OptionModel.MONTE_CARLO,
                                     num_steps=num_steps, num_paths=num_paths, random_seed=random_seed, sens_degree=sens_degree - 1)
            theta = -self.sensitivity('mat', 0.001, model=OptionModel.MONTE_CARLO,
                                      num_steps=num_steps, num_paths=num_paths, random_seed=random_seed, sens_degree=sens_degree - 1)
            rho = self.sensitivity('riskless_rate', 0.0001, model=OptionModel.MONTE_CARLO,
                                   num_steps=num_steps, num_paths=num_paths, random_seed=random_seed, sens_degree=sens_degree - 1)
            vega = self.sensitivity('vol', 0.001, model=OptionModel.MONTE_CARLO,
                                    num_steps=num_steps, num_paths=num_paths, random_seed=random_seed, sens_degree=sens_degree - 1)
        if sens_degree >= 2:
            gamma = self.sensitivity('spot0', self.spot0 * 0.05, opt_measure=OptionMeasure.DELTA, model=OptionModel.MONTE_CARLO,
                                     num_steps=num_steps, num_paths=num_paths, random_seed=random_seed, sens_degree=sens_degree - 1)

        return {
            OptionMeasure.VALUE: val,
            OptionMeasure.DELTA: delta,
            OptionMeasure.THETA: theta,
            OptionMeasure.RHO: rho,
            OptionMeasure.VEGA: vega,
            OptionMeasure.GAMMA: gamma,
        }

    def sensitivity(self, opt_param, opt_param_bump, opt_measure=OptionMeasure.VALUE,
                    model=OptionModel.BINOMIAL_TREE, **model_kwargs):
        up_opt = self.copy()
        setattr(up_opt, opt_param, getattr(up_opt, opt_param) + opt_param_bump)
        down_opt = self.copy()
        setattr(down_opt, opt_param, getattr(down_opt, opt_param) - opt_param_bump)
        up_measure = up_opt.run_model(model=model, **model_kwargs)[opt_measure]
        down_measure = down_opt.run_model(model=model, **model_kwargs)[opt_measure]
        return (up_measure - down_measure) / (2 * opt_param_bump)


def is_number(x):
    return (x is not None) and \
           (isinstance(x, int) or isinstance(x, float) or isinstance(x, numbers.Integral))
