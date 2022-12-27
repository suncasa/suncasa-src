#
# Based on function called in VLA observing script /home/mchost/evla/include/fit_planet_position.py
#
# Original author bjb @ nrao in summer 2012
#
# python functions to fit planet ephemerides and return polynomials
# for ra, dec, and distance.
#
# Slightly modified by Bin Chen @ NJIT on Dec 27, 2022

from math import fabs, sqrt, pow
import numpy as np
#
# this is the singular value decomposition stuff.  it should come from
# numpy, but jython doesn't include it, so i just went ahead and did it
# by hand.  if jython ever includes numpy, then do the following:
#
# from numpy import *
# from scipy import array
#
# as an aside, i really only need searchsorted, polyfit, poly1d, and 
# tolist from numpy, but i can't figure out how to import tolist without 
# importing the whole damn module ><.
#
# then, in the code below, look for the comments related to numpy and
# make the fixes there as well.  also, you don't need the bisect import
# here, because searchsorted does it for you.
#

import bisect


def fpoly(x, np):
    p = [1.0]
    for ii in range(1, np):
        p.append(p[ii - 1] * x)
    return p


def svdfit(xx, yy, sig, ma):
    aa = []
    bb = []
    for ii in range(len(xx)):
        afunc = fpoly(xx[ii], ma)
        tmp = 1.0 / sig[ii]
        taa = []
        for jj in range(ma):
            taa.append(afunc[jj] * tmp)
        bb.append(yy[ii] * tmp)
        aa.append(taa)
    (uu, ww, vv) = svd(aa)
    wmax = max(ww)
    thresh = 2.0e-12 * wmax
    for ii in range(ma):
        if ww[ii] < thresh:
            ww[ii] = 0.0
    aa = svbksb(uu, ww, vv, bb)
    chisq = 0.0
    for ii in range(len(xx)):
        afunc = fpoly(xx[ii], ma)
        model = 0.0
        for jj in range(ma):
            model += aa[jj] * afunc[jj]
        chisq += ((yy[ii] - model) / sig[ii]) * ((yy[ii] - model) / sig[ii])
    return (aa, uu, vv, ww, chisq)


# Almost exact translation of the ALGOL SVD algorithm published in
# Numer. Math. 14, 403-420 (1970) by G. H. Golub and C. Reinsch
#
# Copyright (c) 2005 by Thomas R. Metcalf, helicity314-stitch <at> yahoo <dot> com
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# Pure Python SVD algorithm.
# Input: 2-D list (m by n) with m >= n
# Output: U,W V so that A = U*W*VT
#    Note this program returns V not VT (=transpose(V))
#    On error, a ValueError is raised.
#
# Here is the test case (first example) from Golub and Reinsch
#
# a = [[22.,10., 2.,  3., 7.],
#      [14., 7.,10.,  0., 8.],
#      [-1.,13.,-1.,-11., 3.],
#      [-3.,-2.,13., -2., 4.],
#      [ 9., 8., 1., -2., 4.],
#      [ 9., 1.,-7.,  5.,-1.],
#      [ 2.,-6., 6.,  5., 1.],
#      [ 4., 5., 0., -2., 2.]]
#
# import svd
# import math
# u,w,v = svd.svd(a)
# print w
#
# [35.327043465311384, 1.2982256062667619e-15,
#  19.999999999999996, 19.595917942265423, 0.0]
#
# the correct answer is (the order may vary)
#
# print (math.sqrt(1248.),20.,math.sqrt(384.),0.,0.)
#
# (35.327043465311391, 20.0, 19.595917942265423, 0.0, 0.0)
#
# transpose and matrix multiplication functions are also included
# to facilitate the solution of linear systems.
#
# Version 1.0 2005 May 01

import copy
import math


def svd(a):
    '''Compute the singular value decomposition of array.'''

    # Golub and Reinsch state that eps should not be smaller than the
    # machine precision, ie the smallest number
    # for which 1+e>1.  tol should be beta/e where beta is the smallest
    # positive number representable in the computer.
    eps = 1.e-15  # assumes double precision
    tol = 1.e-64 / eps
    assert 1.0 + eps > 1.0  # if this fails, make eps bigger
    assert tol > 0.0  # if this fails, make tol bigger
    itmax = 50
    u = copy.deepcopy(a)
    m = len(a)
    n = len(a[0])
    #   if __debug__: print 'a is ',m,' by ',n

    if m < n:
        if __debug__: print('Error: m is less than n')
        raise ValueError('SVD Error: m is less than n.')

    e = [0.0] * n  # allocate arrays
    q = [0.0] * n
    v = []
    for k in range(n): v.append([0.0] * n)

    # Householder's reduction to bidiagonal form

    g = 0.0
    x = 0.0

    for i in range(n):
        e[i] = g
        s = 0.0
        l = i + 1
        for j in range(i, m): s += (u[j][i] * u[j][i])
        if s <= tol:
            g = 0.0
        else:
            f = u[i][i]
            if f < 0.0:
                g = math.sqrt(s)
            else:
                g = -math.sqrt(s)
            h = f * g - s
            u[i][i] = f - g
            for j in range(l, n):
                s = 0.0
                for k in range(i, m): s += u[k][i] * u[k][j]
                f = s / h
                for k in range(i, m): u[k][j] = u[k][j] + f * u[k][i]
        q[i] = g
        s = 0.0
        for j in range(l, n): s = s + u[i][j] * u[i][j]
        if s <= tol:
            g = 0.0
        else:
            f = u[i][i + 1]
            if f < 0.0:
                g = math.sqrt(s)
            else:
                g = -math.sqrt(s)
            h = f * g - s
            u[i][i + 1] = f - g
            for j in range(l, n): e[j] = u[i][j] / h
            for j in range(l, m):
                s = 0.0
                for k in range(l, n): s = s + (u[j][k] * u[i][k])
                for k in range(l, n): u[j][k] = u[j][k] + (s * e[k])
        y = abs(q[i]) + abs(e[i])
        if y > x: x = y
    # accumulation of right hand gtransformations
    for i in range(n - 1, -1, -1):
        if g != 0.0:
            h = g * u[i][i + 1]
            for j in range(l, n): v[j][i] = u[i][j] / h
            for j in range(l, n):
                s = 0.0
                for k in range(l, n): s += (u[i][k] * v[k][j])
                for k in range(l, n): v[k][j] += (s * v[k][i])
        for j in range(l, n):
            v[i][j] = 0.0
            v[j][i] = 0.0
        v[i][i] = 1.0
        g = e[i]
        l = i
    # accumulation of left hand transformations
    for i in range(n - 1, -1, -1):
        l = i + 1
        g = q[i]
        for j in range(l, n): u[i][j] = 0.0
        if g != 0.0:
            h = u[i][i] * g
            for j in range(l, n):
                s = 0.0
                for k in range(l, m): s += (u[k][i] * u[k][j])
                f = s / h
                for k in range(i, m): u[k][j] += (f * u[k][i])
            for j in range(i, m): u[j][i] = u[j][i] / g
        else:
            for j in range(i, m): u[j][i] = 0.0
        u[i][i] += 1.0
    # diagonalization of the bidiagonal form
    eps = eps * x
    for k in range(n - 1, -1, -1):
        for iteration in range(itmax):
            # test f splitting
            for l in range(k, -1, -1):
                goto_test_f_convergence = False
                if abs(e[l]) <= eps:
                    # goto test f convergence
                    goto_test_f_convergence = True
                    break  # break out of l loop
                if abs(q[l - 1]) <= eps:
                    # goto cancellation
                    break  # break out of l loop
            if not goto_test_f_convergence:
                # cancellation of e[l] if l>0
                c = 0.0
                s = 1.0
                l1 = l - 1
                for i in range(l, k + 1):
                    f = s * e[i]
                    e[i] = c * e[i]
                    if abs(f) <= eps:
                        # goto test f convergence
                        break
                    g = q[i]
                    h = pythag(f, g)
                    q[i] = h
                    c = g / h
                    s = -f / h
                    for j in range(m):
                        y = u[j][l1]
                        z = u[j][i]
                        u[j][l1] = y * c + z * s
                        u[j][i] = -y * s + z * c
            # test f convergence
            z = q[k]
            if l == k:
                # convergence
                if z < 0.0:
                    # q[k] is made non-negative
                    q[k] = -z
                    for j in range(n):
                        v[j][k] = -v[j][k]
                break  # break out of iteration loop and move on to next k value
            if iteration >= itmax - 1:
                if __debug__: print('Error: no convergence.')
                # should this move on the the next k or exit with error??
                # raise ValueError,'SVD Error: No convergence.'  # exit the program with error
                break  # break out of iteration loop and move on to next k
            # shift from bottom 2x2 minor
            x = q[l]
            y = q[k - 1]
            g = e[k - 1]
            h = e[k]
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y)
            g = pythag(f, 1.0)
            if f < 0:
                f = ((x - z) * (x + z) + h * (y / (f - g) - h)) / x
            else:
                f = ((x - z) * (x + z) + h * (y / (f + g) - h)) / x
            # next QR transformation
            c = 1.0
            s = 1.0
            for i in range(l + 1, k + 1):
                g = e[i]
                y = q[i]
                h = s * g
                g = c * g
                z = pythag(f, h)
                e[i - 1] = z
                c = f / z
                s = h / z
                f = x * c + g * s
                g = -x * s + g * c
                h = y * s
                y = y * c
                for j in range(n):
                    x = v[j][i - 1]
                    z = v[j][i]
                    v[j][i - 1] = x * c + z * s
                    v[j][i] = -x * s + z * c
                z = pythag(f, h)
                q[i - 1] = z
                c = f / z
                s = h / z
                f = c * g + s * y
                x = -s * g + c * y
                for j in range(m):
                    y = u[j][i - 1]
                    z = u[j][i]
                    u[j][i - 1] = y * c + z * s
                    u[j][i] = -y * s + z * c
            e[l] = 0.0
            e[k] = f
            q[k] = x
            # goto test f splitting

    # vt = transpose(v)
    # return (u,q,vt)
    return (u, q, v)


def pythag(a, b):
    absa = abs(a)
    absb = abs(b)
    if absa > absb:
        return absa * math.sqrt(1.0 + (absb / absa) ** 2)
    else:
        if absb == 0.0:
            return 0.0
        else:
            return absb * math.sqrt(1.0 + (absa / absb) ** 2)


def svbksb(uu, ww, vv, bb):
    ma = len(uu[0])
    nn = len(uu)
    tmp = []
    for ii in range(ma):
        ss = 0.0
        if ww[ii] != 0.0:
            for jj in range(nn):
                ss += uu[jj][ii] * bb[jj]
            ss /= ww[ii]
        tmp.append(ss)
    xx = []
    for ii in range(ma):
        ss = 0.0
        for jj in range(ma):
            ss += vv[ii][jj] * tmp[jj]
        xx.append(ss)
    return xx


def svdvar(vv, ww, ma):
    wti = []
    for ii in range(ma):
        wti.append(0.0)
        if ww[ii] != 0.0:
            wti[ii] = 1.0 / (ww[ii] * ww[ii])
    cvm = [[0 for col in range(ma)] for row in range(ma)]
    for ii in range(ma):
        for jj in range(ii + 1):
            sum = 0.0
            for kk in range(ma):
                sum += vv[ii][kk] * vv[jj][kk] * wti[kk]
            cvm[ii][jj] = sum
            cvm[jj][ii] = sum
    return cvm


def fit_planet_positions(times, ras, decs, start_time=None, end_time=None, distances=None, allowed_error=0.01):
    '''
    find a fitting polynomial for an ephemeris table.

    inputs:
        times = list of times at which the ephemeris positions are
                tabulated.  assumed sorted in ascending order.
                (MJD).
        ras = list of right ascensions at those times (radians).
        decs = list of declinations at those times (radians).
    Optional inputs:
        start_time = the start time of the SB (MJD).
        end_time = the end time of the SB (MJD).
        distances = list of distances at those times (AU).
        allowed_error = the allowed error in the fitting polynomials
                        for ra and dec from the tabulated values
                        (asec).

    returned is a list, first element is the return status:
        0 -> success
        1-7 -> Warning: did not converge to required accuracy.  
               1 -> right ascension only didn't converge
               2 -> declination only didn't converge
               3 -> right ascension and declination didn't converge
               4 -> distance only didn't converge
               5 -> right ascension and distance didn't converge
               6 -> declination and distance didn't converge
               7 -> none of the three converged
               (note that in this case the best fitted polynomials 
               are still returned [N.B. i should return the error 
               somewhere...].)
        8 -> Error: the time range from start_time to end_time is
             not completely contained in the tabulated times.
    second element is the time t0 to which the polynomial is
        referenced.
    third element is the list of right ascension coefficients.
    fourth element is the list of declination coefficients.
    fifth element is the list of distance coefficients.

    bjb
    nrao
    summer 2012
    '''

    NPOLYMAX = 7
    ASEC2RAD = 2.0626480624710e5

    times = np.array(times)
    ras = np.array(ras)
    decs = np.array(decs)

    #
    # check the times
    #
    if start_time or end_time:
        if start_time < times[0] or end_time > times[-1]:
            result = [8]
            return result

    result = [0]

    #
    # convert from arcseconds to radians
    #
    allowed_error /= ASEC2RAD

    #
    # carve out the sub-lists for the requested time range
    #
    # if numpy ever gets into jython, use the following
    # instead of what is below...
    #
    #   start_index = searchsorted(times, start_time)
    #   end_index = searchsorted(times, end_time)
    #
    # use the bisect functions in lieu of numpy/searchsorted
    #
    if start_time:
        start_index = bisect.bisect_left(times, start_time)
        start_index = max(start_index - 3, 0)
    else:
        start_index = 0
    if end_time:
        end_index = bisect.bisect_right(times, end_time)
    else:
        end_index = len(times) - 1

    end_index = min(end_index + 3, len(times) - 1)
    times_slice = times[start_index:end_index]
    ras_slice = ras[start_index:end_index]
    decs_slice = decs[start_index:end_index]

    if distances:
        distances_slice = distances[start_index:end_index]

    #
    # reference times to the center time, otherwise the
    # polynomials have huge coefficients and potentially
    # blow up
    #
    t0 = times_slice[int(len(times_slice) / 2)]
    #print('length of the time array used for polyfit', len(times_slice))
    result.append(t0)
    times_slice -= t0
    #for ii in range(len(times_slice)):
    #    times_slice[ii] -= t0

    #
    # make the corresponding scipy/numpy array quantities
    #
    #   time_array = array(times_slice)
    #   ra_array = array(ras_slice)
    #   dec_array = array(decs_slice)
    #   distance_array = array(distances_slice)
    time_array = times_slice
    ra_array = ras_slice
    dec_array = decs_slice
    if distances:
        distance_array = distances_slice

    #
    # right ascension fit
    #
    #
    # for non-numpy version, need errors.  just take them all equal.
    #
    sig = []
    for ii in range(len(ra_array)):
        sig.append(ra_array[0] / 1000.0)
    max_ra_difn = allowed_error + 1.0
    npoly_ra = 2
    while max_ra_difn > allowed_error and npoly_ra <= NPOLYMAX:
        #
        # if numpy ever gets into jython...
        #
        #       a_ra = polyfit (time_array, ra_array, npoly_ra)
        #       ra_fit = poly1d(a_ra)
        #
        # polyfit returns the coefficients in reverse order from normal,
        # so reverse them back...
        #
        #       a_ra = a_ra[::-1]
        #       max_ra_difn = 0.0
        #       for ii in range(len(time_array)):
        #           model_ra = ra_fit(time_array[ii])
        #           max_ra_difn = max (max_ra_difn, fabs(model_ra-ra_array[ii]))
        #
        # this is the numpy replacement, done by hand...
        #
        (a_ra, uu, vt, ww, chisq) = svdfit(time_array, ra_array, sig, npoly_ra)
        #
        # i never use the covariance matrix, so could comment this out,
        # but i leave it in, in case i ever want to check on the significance
        # of the coefficients (compare them against their error) to see whether
        # there is a problem.
        #
        cvm = svdvar(vt, ww, npoly_ra)
        max_ra_difn = 0.0
        for ii in range(len(time_array)):
            model_ra = 0.0
            for jj in range(npoly_ra):
                model_ra += a_ra[jj] * pow(time_array[ii], jj)
            max_ra_difn = max(max_ra_difn, fabs(model_ra - ra_array[ii]))
        npoly_ra += 1
    npoly_ra -= 1
    if npoly_ra == NPOLYMAX and max_ra_difn > allowed_error:
        result[0] = 1
    #
    # the numpy version - have to convert from numpy array to list
    #
    #   result.append(a_ra.tolist())
    result.append(a_ra)

    #
    # declination fit
    #
    sig = []
    for ii in range(len(dec_array)):
        sig.append(dec_array[0] / 1000.0)
    max_dec_difn = allowed_error + 1.0
    npoly_dec = 1
    while max_dec_difn > allowed_error and npoly_dec <= NPOLYMAX:
        # numpy
        #       a_dec = polyfit (time_array, dec_array, npoly_dec)
        #       dec_fit = poly1d(a_dec)
        #       a_dec = a_dec[::-1]
        #       max_dec_difn = 0.0
        #       for ii in range(len(time_array)):
        #           model_dec = dec_fit(time_array[ii])
        #           max_dec_difn = max (max_dec_difn, fabs(model_dec-dec_array[ii]))
        (a_dec, uu, vt, ww, chisq) = svdfit(time_array, dec_array, sig, npoly_dec)
        cvm = svdvar(vt, ww, npoly_dec)
        max_dec_difn = 0.0
        for ii in range(len(time_array)):
            model_dec = 0.0
            for jj in range(npoly_dec):
                model_dec += a_dec[jj] * pow(time_array[ii], jj)
            max_dec_difn = max(max_dec_difn, fabs(model_dec - dec_array[ii]))
        npoly_dec += 1
    npoly_dec -= 1
    if npoly_dec == NPOLYMAX and max_dec_difn > allowed_error:
        result[0] += 2
    #   result.append(a_dec.tolist())
    result.append(a_dec)

    #
    # distance fit.  for distance, set the allowed error to ~1.5 km...
    #
    if distances:
        sig = []
        for ii in range(len(distance_array)):
            sig.append(distance_array[0] / 1000.0)
        allowed_error = 1.0e-8
        max_distance_difn = allowed_error + 1.0
        npoly_distance = 1
        while max_distance_difn > allowed_error and npoly_distance <= NPOLYMAX:
            # numpy
            #       a_distance = polyfit (time_array, distance_array, npoly_distance)
            #       distance_fit = poly1d(a_distance)
            #       a_distance = a_distance[::-1]
            #       max_distance_difn = 0.0
            #       for ii in range(len(time_array)):
            #           model_distance = distance_fit(time_array[ii])
            #           max_distance_difn = max (max_distance_difn, fabs(model_distance-distance_array[ii]))
            (a_distance, uu, vt, ww, chisq) = svdfit(time_array, distance_array, sig, npoly_distance)
            cvm = svdvar(vt, ww, npoly_distance)
            max_distance_difn = 0.0
            for ii in range(len(time_array)):
                model_distance = 0.0
                for jj in range(npoly_distance):
                    model_distance += a_distance[jj] * pow(time_array[ii], jj)
                max_distance_difn = max(max_distance_difn, fabs(model_distance - distance_array[ii]))
            npoly_distance += 1
        npoly_distance -= 1
        if npoly_distance == NPOLYMAX and max_distance_difn > allowed_error:
            result[0] += 4
        #   result.append(a_distance.tolist())
        result.append(a_distance)

    return result
