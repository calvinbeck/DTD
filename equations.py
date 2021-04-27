#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Equations:
Betaversion: Will be improved in near future
"""
import numpy as np


def d2Tdx2(T_xm, T_x, T_xp, Delta_x):
    h_x = Delta_x[2] - Delta_x[1]
    h_xm = Delta_x[1] - Delta_x[0]
    numerator = 2 * (T_xp + h_x / h_xm * T_xm - (1 + h_x / h_xm) * T_x)
    denominator = h_x * h_xm * (1 + h_x / h_xm)
    return numerator/denominator


def dTdt(T_tm, T_tp, Delta_t):
    return (T_tp-T_tm)/(2*Delta_t)


def d2Tdx2_error(T_xm, T_x, T_xp, Delta_x, T_error, Delta_x_error):
    h_x = Delta_x[2] - Delta_x[1]
    h_xm = Delta_x[1] - Delta_x[0]
    Delta_x_error = np.sqrt(2)*Delta_x_error
    sigma_T = T_error**2*(1/(h_xm**2*(h_x+h_xm)**2) + 1
                          / (h_x**2*(h_x+h_xm)**2) + 1/(h_x**2*h_xm**2))
    sigma_x1 = Delta_x_error**2*(T_xm*h_x**2 + 2*T_xm*h_x*h_xm - T_x*h_x**2
                                 - 2*T_x*h_x*h_xm-T_x*h_xm**2
                                 + T_xp*h_xm**2)**2 / (h_x**2*h_xm**4
                                 * (h_x+h_xm)**4)
    sigma_x2 = Delta_x_error**2*(T_xm*h_x**2 - T_x*h_x**2 - 2*T_x*h_x*h_xm
                                 - T_x*h_xm**2 + 2*T_xp*h_x*h_xm
                                 + T_xp*h_xm**2)**2/(h_x**4*h_xm**2
                                 * (h_x+h_xm)**4)
    sigma_d2Tdx2 = 2 * np.sqrt(sigma_T + sigma_x1 + sigma_x2)
    return sigma_d2Tdx2


def dTdt_error(T_tm, T_tp, Delta_t, T_error):
    sigma_dTdt = T_error/(np.sqrt(2)*Delta_t)
    return sigma_dTdt


"""
print(d2Tdx2(215, 214, 213.1, [0.,0.11,0.21]))  # 0.8658008657946646
print(d2Tdx2_error(215, 214, 213.1, [0.,0.11,0.21], 0.1, 0.01)) # 27.694318768

print(dTdt(215, 213.1, 60)) # -0.01583333333333338
print(dTdt_error(215, 213.1, 60, 0.1))  # 0.00117851130198

"""