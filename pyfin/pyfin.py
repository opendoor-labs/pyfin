# -*- coding: utf-8 -*-
import abc
from math import exp
import math

import scipy.optimize
import scipy.stats


def enum(**enums):
    return type('Enum', (), enums)

OptionType = enum(CALL='call', PUT='put')
OptionExerciseType = enum(EUROPEAN='european', AMERICAN='american')
OptionModel = enum(BLACK_SCHOLES='black_scholes', BINOMIAL_TREE='binomial_tree')
OptionMeasure = enum(VALUE='value', DELTA='delta', THETA='theta', RHO='rho', VEGA='vega', GAMMA='gamma')

DEFAULT_BINOMIAL_TREE_NUM_STEPS = 25


class Instrument(object):
    @abc.abstractmethod
    def run_model(self):
        """Calculate a measures (i.e. theoretical value & greeks) for this instrument"""


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
    def imply_volatility(premium, *args, **kwargs):
        def obj_fn(vol_guess):
            kwargs['vol'] = vol_guess
            val = Option(*args, **kwargs).run_model()[OptionMeasure.VALUE]
            return val - premium
        try:
            return scipy.optimize.bisect(obj_fn, 0.01, 0.99, xtol=0.0025)
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
        else:  # model == OptionModel.BINOMIAL_MODEL:
            result = self.binomial_tree(**kwargs)

        self.model_cache[model] = result
        self.model_cache_param_hashes[model] = curr_param_hash
        return result

    def black_scholes(self):
        assert self.exer_type == OptionExerciseType.EUROPEAN, "Black-Scholes does not support early exercise"

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

    def binomial_tree(self, num_steps=None, is_greeks=True):
        val_num_steps = num_steps or DEFAULT_BINOMIAL_TREE_NUM_STEPS
        val_cache = {}

        dt = float(self.mat) / val_num_steps
        up_fact = exp(self.vol * (dt ** 0.5))
        down_fact = 1.0 / up_fact
        prob_up = (exp((self.riskless_rate - self.yield_) * dt) - down_fact) / (up_fact - down_fact)

        def spot_price(num_ups, num_downs):
            return self.spot0 * (up_fact ** (num_ups - num_downs))

        def value_helper(n, num_ups, num_downs):
            value_cache_key = (n, num_ups, num_downs)
            if value_cache_key not in val_cache:
                spot = spot_price(num_ups, num_downs)
                if self.opt_type == OptionType.CALL:
                    exer_profit = max(0, spot - self.strike)
                else:  # opt_type == OptionType.PUT:
                    exer_profit = max(0, self.strike - spot)

                if n >= val_num_steps:
                    val = exer_profit
                else:
                    fv = prob_up * value_helper(n + 1, num_ups + 1, num_downs) + \
                        (1 - prob_up) * value_helper(n + 1, num_ups, num_downs + 1)
                    pv = exp(-self.riskless_rate * dt) * fv
                    if self.exer_type == OptionExerciseType.AMERICAN:
                        val = max(pv, exer_profit)
                    else:  # exer_type == OptionExerciseType.EUROPEAN
                        val = pv

                val_cache[value_cache_key] = val
            return val_cache[value_cache_key]

        val = value_helper(0, 0, 0)
        delta, theta, rho, vega, gamma = None, None, None, None, None
        if is_greeks:
            delta = (value_helper(1, 1, 0) - value_helper(1, 0, 1)) / (self.spot0 * up_fact - self.spot0 * down_fact)
            theta = (value_helper(2, 1, 1) - val) / (2 * dt)

            rho = self.value_sensitivity('riskless_rate', 0.0025, model=OptionModel.BINOMIAL_TREE)
            vega = self.value_sensitivity('vol', 0.05, model=OptionModel.BINOMIAL_TREE)

            delta_up = (value_helper(2, 2, 0) - value_helper(2, 1, 1)) / (self.spot0 * up_fact - self.spot0 * down_fact)
            delta_down = (value_helper(2, 1, 1) - value_helper(2, 0, 2)) / (self.spot0 * up_fact - self.spot0 * down_fact)
            gamma = (delta_up - delta_down) / ((self.spot0 * up_fact * up_fact - self.spot0 * down_fact * down_fact) / 2)

        return {
            OptionMeasure.VALUE: val,
            OptionMeasure.DELTA: delta,
            OptionMeasure.THETA: theta,
            OptionMeasure.RHO: rho,
            OptionMeasure.VEGA: vega,
            OptionMeasure.GAMMA: gamma,
        }

    def value_sensitivity(self, opt_param, bump, model=OptionModel.BINOMIAL_TREE):
        up_opt = self.copy()
        setattr(up_opt, opt_param, getattr(up_opt, opt_param) + bump)
        down_opt = self.copy()
        setattr(down_opt, opt_param, getattr(down_opt, opt_param) - bump)
        up_val = up_opt.run_model(model=model, is_greeks=False)[OptionMeasure.VALUE]
        down_val = down_opt.run_model(model=model, is_greeks=False)[OptionMeasure.VALUE]
        return (up_val - down_val) / (2 * bump)
