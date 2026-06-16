# SPDX-License-Identifier: MIT

"""
Wind Turbine Annual Energy Production Calculator
================================================

Author:
    Janusz Telega

Email:
    janusz.telega@imp.gda.pl

Affiliation:
    Institute of Fluid-Flow Machinery,
    Polish Academy of Sciences,
    Gda?sk, Poland

Date:
    2026-06-15

Version:
    1.0

Description:
    This script estimates the monthly and annual electrical energy production
    of a small wind turbine using:

        1. a turbine power curve P(v),
        2. monthly Weibull wind speed distributions,
        3. discrete numerical integration over wind speed bins.

    The turbine power curve should be provided as electrical output power
    versus wind speed. It may be obtained experimentally, for example from
    wind tunnel measurements.

Main assumptions:
    - Wind speed is described by monthly Weibull distributions.
    - The Weibull shape parameter is denoted as k [-].
    - The Weibull scale parameter is denoted as c [m/s].
    - The turbine power curve is interpolated linearly between given points.
    - Power outside the defined operating range is assumed to be zero when
      np.interp(..., left=0.0, right=0.0) is used.
    - Air density correction is not included in this basic version.

Units:
    - Wind speed:        m/s
    - Turbine power:     W
    - Mean power:        W
    - Energy:            Wh or kWh
    - Time:              h

Input data:
    k_month:
        Monthly Weibull shape parameters [-].

    c_month:
        Monthly Weibull scale parameters [m/s].

    v_curve:
        Wind speed points of the turbine power curve [m/s].

    p_curve_W:
        Electrical power output of the turbine at v_curve points [W].

Output:
    - Mean monthly power [W]
    - Monthly energy production [Wh, kWh]
    - Annual energy production [Wh, kWh]

Reference wind data:
    Monthly Weibull parameters may be replaced by local wind statistics.

    Example reference dataset:
    W. Zhang, J. Harff, R. Schneider,
    "Analysis of 50-year wind data of the southern Baltic Sea for modelling
    coastal morphological evolution ? a case study from the Darss-Zingst
    Peninsula."

License:
    MIT License

    Copyright (c) 2026 Janusz Telega

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
    IN THE SOFTWARE.
"""


import numpy as np

# Monthly Weibull parameters
# After Table 2 in: W. Zhang, J. Harff, R. Schneider, Analysis of 50-year wind data of the southern
# Baltic Sea for modelling coastal morphological evolution ? a case study from 
# the Darss-Zingst Peninsula, OCEANOLOGIA, 53 (1-TI), 2011. pp. 489?518.


k_month = np.array([2.44, 2.41, 2.29, 2.24, 2.22, 2.17, 2.12, 2.20, 2.27, 2.40, 2.51, 2.52])

c_month = np.array([9.40, 8.97, 8.52, 7.16, 6.70, 6.66, 6.84, 6.87, 7.61, 8.89, 9.29, 9.37])

hours_month = np.array([
    31*24, 28*24, 31*24, 30*24,
    31*24, 30*24, 31*24, 31*24,
    30*24, 31*24, 30*24, 31*24
])

# Turnine characteristic curve : velocity [m/s], power [W]
# v_curve = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 25])
# p_curve_W = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1])  # W
v_curve = np.array([3, 4, 5, 6, 7, 8, 9, 10, 11, 12])
p_curve_W = np.array([1, 1, 1, 1, 1, 1, 1, 1, 1, 1])  # W

def weibull_cdf(v, k, c):
    return 1.0 - np.exp(-(v / c)**k)

dv = 0.1
v_edges = np.arange(0, 40 + dv, dv)
v_mid = 0.5 * (v_edges[:-1] + v_edges[1:])

energy_month_Wh = []
mean_power_month_W = []

for k, c, hours in zip(k_month, c_month, hours_month):
    prob = weibull_cdf(v_edges[1:], k, c) - weibull_cdf(v_edges[:-1], k, c)

    # Important: right=0 means that the last point of characteristic is a cut-off value; for highger values of velocity, zero power is generated (the Weibull goes to +infty).
    # If there is suppression for constant/rated power, the 'right=0.0' should be erased
    power_mid_W = np.interp(v_mid, v_curve, p_curve_W, left=0.0, right=0.0)

    mean_power_W = np.sum(power_mid_W * prob)   # W
    energy_Wh = mean_power_W * hours            # Wh

    mean_power_month_W.append(mean_power_W)
    energy_month_Wh.append(energy_Wh)

energy_month_Wh = np.array(energy_month_Wh)
mean_power_month_W = np.array(mean_power_month_W)

energy_year_Wh = np.sum(energy_month_Wh)
energy_year_kWh = energy_year_Wh / 1000.0

print(f"Annual production: {energy_year_Wh:.0f} Wh")
print(f"Annual production: {energy_year_kWh:.2f} kWh")

for i, (p_mean, e_Wh) in enumerate(zip(mean_power_month_W, energy_month_Wh), start=1):
    print(f"Month {i:2d}: mean power = {p_mean:.2f} W, energy = {e_Wh:.1f} Wh = {e_Wh/1000:.3f} kWh")