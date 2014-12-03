# pyfin
-------

[<img src="https://badge.fury.io/py/pyfin.png">](http://badge.fury.io/py/pyfin)
[<img src="https://travis-ci.org/opendoor-labs/pyfin.png?branch=master">](https://travis-ci.org/opendoor-labs/pyfin)
[<img src="https://pypip.in/d/pyfin/badge.png">](https://pypi.python.org/pypi/pyfin)

###### Basic options pricing in Python

&ldquo;Oh cool. Probably a little easier than spinning up the QuantLib stack.&rdquo; &mdash; [Wes McKinney](https://github.com/wesm), creator of [Pandas](https://github.com/pydata/pandas)


### Features
-------

* Option valuation w/ Black-Scholes, lattice (binomial tree), and Monte Carlo simulation models.
* Basic Greeks calculation (delta, theta, rho, vega, gamma) across each valuation model.
* Discrete dividends support in the lattice (binomial tree) and Monte Carlo simulation models.
* Early exercise (American options) support in Monte Carlo simulation through the Longstaff-Schwartz technique.
* Minimal dependencies, just Numpy & SciPy.
* Free software, released under the MIT license.
