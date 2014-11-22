===============================
pyfin
===============================

.. image:: https://badge.fury.io/py/pyfin.png
    :target: http://badge.fury.io/py/pyfin

.. image:: https://travis-ci.org/opendoor-labs/pyfin.png?branch=master
        :target: https://travis-ci.org/opendoor-labs/pyfin

.. image:: https://pypip.in/d/pyfin/badge.png
        :target: https://pypi.python.org/pypi/pyfin


Basic options pricing in Python

* Free software: MIT license

Features
--------

* Option valuation w/ Black-Scholes, lattice (binomial tree), and Monte Carlo simulation models.
* Basic Greeks calculation (delta, theta, rho, vega, gamma) across each valuation model.
* Discrete dividends support in the lattice (binomial tree) and Monte Carlo simulation models.
* Early exercise (American options) support in Monte Carlo simulation through the Longstaff-Schwartz technique.
* Minimal dependencies, just Numpy & SciPy.
