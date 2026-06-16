% SPDX-License-Identifier: MIT
%
% Wind Turbine Annual Energy Production Calculator
% =================================================
%
% Author:
%     Janusz Telega
%
% Email:
%     janusz.telega@imp.gda.pl
%
% Affiliation:
%     Institute of Fluid-Flow Machinery,
%     Polish Academy of Sciences,
%     Gdańsk, Poland
%
% Date:
%     2026-06-15
%
% Version:
%     1.0
%
% Description:
%     This script estimates the monthly and annual electrical energy
%     production of a small wind turbine using:
%
%         1. a turbine power curve P(v),
%         2. monthly Weibull wind speed distributions,
%         3. discrete numerical integration over wind speed bins.
%
%     The turbine power curve should be provided as electrical output power
%     versus wind speed. It may be obtained experimentally, for example from
%     wind tunnel measurements.
%
% Main assumptions:
%     - Wind speed is described by monthly Weibull distributions.
%     - The Weibull shape parameter is denoted as k [-].
%     - The Weibull scale parameter is denoted as c [m/s].
%     - The turbine power curve is interpolated linearly between given points.
%     - Power outside the defined operating range is assumed to be zero.
%     - Air density correction is not included in this basic version.
%
% Units:
%     - Wind speed:        m/s
%     - Turbine power:     W
%     - Mean power:        W
%     - Energy:            Wh or kWh
%     - Time:              h
%
% License:
%     MIT License
%
%     Copyright (c) 2026 Janusz Telega
%
% -------------------------------------------------------------------------

clear;
clc;

%% Monthly Weibull parameters
% After Table 2 in: W. Zhang, J. Harff, R. Schneider, Analysis of 50-year wind data of the southern
% Baltic Sea for modelling coastal morphological evolution – a case study from 
% the Darss-Zingst Peninsula, OCEANOLOGIA, 53 (1-TI), 2011. pp. 489–518.

k_month = [2.44, 2.41, 2.29, 2.24, 2.22, 2.17, ...
           2.12, 2.20, 2.27, 2.40, 2.51, 2.52];

c_month = [9.40, 8.97, 8.52, 7.16, 6.70, 6.66, ...
           6.84, 6.87, 7.61, 8.89, 9.29, 9.37];

hours_month = [31*24, 28*24, 31*24, 30*24, ...
               31*24, 30*24, 31*24, 31*24, ...
               30*24, 31*24, 30*24, 31*24];

%% Turbine power curve
% Wind speed [m/s]

v_curve = [3, 4, 5, 6, 7, 8, 9, 10, 11, 12];

% Electrical power output [W]

p_curve_W = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1];

%% Wind speed discretization

dv = 0.1;

v_edges = 0:dv:40;
v_mid = 0.5 * (v_edges(1:end-1) + v_edges(2:end));

%% Interpolated turbine power
%
% MATLAB interp1 with extrapolation value = 0 reproduces:
% np.interp(..., left=0.0, right=0.0)

power_mid_W = interp1( ...
    v_curve, ...
    p_curve_W, ...
    v_mid, ...
    'linear', ...
    0);

%% Preallocation

energy_month_Wh = zeros(1,12);
mean_power_month_W = zeros(1,12);

%% Monthly calculations

for m = 1:12

    k = k_month(m);
    c = c_month(m);
    hours = hours_month(m);

    % Probability in each velocity bin
    prob = weibull_cdf(v_edges(2:end), k, c) - ...
           weibull_cdf(v_edges(1:end-1), k, c);

    % Mean power [W]
    mean_power_W = sum(power_mid_W .* prob);

    % Monthly energy [Wh]
    energy_Wh = mean_power_W * hours;

    mean_power_month_W(m) = mean_power_W;
    energy_month_Wh(m) = energy_Wh;

end
% figure
% plot(v_mid,prob)
% xlabel( 'V [m/s]')
% title('Probability distribution function')
%% Annual production

energy_year_Wh = sum(energy_month_Wh);
energy_year_kWh = energy_year_Wh / 1000;

fprintf('Annual production: %.0f Wh\n', energy_year_Wh);
fprintf('Annual production: %.2f kWh\n\n', energy_year_kWh);

%% Monthly results

for m = 1:12

    fprintf(['Month %2d: mean power = %.2f W, ' ...
             'energy = %.1f Wh = %.3f kWh\n'], ...
             m, ...
             mean_power_month_W(m), ...
             energy_month_Wh(m), ...
             energy_month_Wh(m)/1000);

end

%% Local function

function F = weibull_cdf(v, k, c)
%WEIBULL_CDF Weibull cumulative distribution function
%
% F(v) = 1 - exp(-(v/c)^k)

    F = 1 - exp(-(v ./ c).^k);

end