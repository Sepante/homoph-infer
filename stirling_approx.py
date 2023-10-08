#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 13 15:13:50 2021

@author: sepante
"""

import numpy as np
from decimal import Decimal
#from scipy.special import factorial, comb

def factorial( n ):
    if n == 0:
        return Decimal(1)
    else:
        return np.sqrt( Decimal(2) * Decimal(np.pi) * n ) * ((n / Decimal(np.e) ) ** n)

def comb( n, k ):
    if (n == 0) or (k == 0):
        return 1
    elif k > n:
        return 0
    else:
        return factorial( n ) / factorial( k ) / factorial( n - k )


def total_pairing_num(total_stubs):
    if (total_stubs % 2 != 0):
        print('number of stubs should be even.')
    else:
        numin = factorial( total_stubs)
        denom_1 = factorial(total_stubs//2)
        denom_2 = Decimal(2 ** (total_stubs // 2))


        return numin / denom_1 / denom_2

def pairing_cond_num_asym( stub_kind_1, stub_kind_2, phi, h ): #phi: number of heterogeneous pairs
    stub_kind_1 = Decimal( stub_kind_1 )
    stub_kind_2 = Decimal( stub_kind_2 )

    pairing_num = total_pairing_num( stub_kind_1 - phi) \
 * total_pairing_num( stub_kind_2 - phi)\
 * comb( stub_kind_1, phi)\
 * comb( stub_kind_2, phi)\
 * factorial(phi)
    # print(h)
    h = Decimal(h)
    # print(h)
    pairing_prob = pairing_num * ( (1-h) / h ) ** phi
    return pairing_prob

def pairing_cond_num_asym_norm_fact( stub_kind_1, stub_kind_2, h ): #phi: number of heterogeneous pairs
    s_minority = np.min([stub_kind_1, stub_kind_2])

    min_phi = stub_kind_1 % 2 # if stub groups are odd, we need at least 1 pairing edge, otherwise 0 is the minimum

    phi_seq = range(min_phi, s_minority, 2)

    norm_factor = 0

    for phi in phi_seq:
        pairing_prob = pairing_cond_num_asym( stub_kind_1, stub_kind_2, phi, h )
#         print(pairing_prob)
        norm_factor += pairing_prob


    return norm_factor

def pairing_cond_num_asym_asym_h( stub_kind_1, stub_kind_2, phi, h_11, h_22 ): #phi: number of heterogeneous pairs
    h_11 = Decimal(h_11)
    h_22 = Decimal(h_22)
    inter_community_edges_configs = 0
    mue_seq = [Decimal(mue) for mue in range(phi)]
    for mue in mue_seq:
        inter_community_edges_configs += comb(phi, mue) * (1-h_11) ** (mue) * (1-h_22) ** (phi - mue)



    stub_kind_1 = Decimal( stub_kind_1 )
    stub_kind_2 = Decimal( stub_kind_2 )

    pairing_prob =\
   total_pairing_num( stub_kind_1 - phi)\
 * total_pairing_num( stub_kind_2 - phi)\
 * comb( stub_kind_1, phi)\
 * comb( stub_kind_2, phi)\
 * factorial(phi)\
 * (2*h_11) ** ((stub_kind_1 - phi)/2)\
 * (2*h_22) ** ((stub_kind_2 - phi)/2)\
 * inter_community_edges_configs


#     h = Decimal(h)
#     pairing_prob = pairing_num * ( (1-h) / h ) ** phi
    return pairing_prob

def pairing_cond_num_asym_asym_norm_fact( stub_kind_1, stub_kind_2, h_11, h_22 ): #phi: number of heterogeneous pairs
    s_minority = np.min([stub_kind_1, stub_kind_2])

    min_phi = stub_kind_1 % 2 # if stub groups are odd, we need at least 1 pairing edge, otherwise 0 is the minimum

    phi_seq = range(min_phi, s_minority, 2)

    norm_factor = 0

    for phi in phi_seq:
        pairing_prob = pairing_cond_num_asym_asym_h( stub_kind_1, stub_kind_2, phi, h_11, h_22 )
#         print(pairing_prob)
        norm_factor += pairing_prob


    return norm_factor



def pairing_cond_num_asym_asym_h_mue_conditioned( stub_kind_1, stub_kind_2, phi, h_11, h_22, mue ): #phi: number of heterogeneous pairs
    h_11 = Decimal(h_11)
    h_22 = Decimal(h_22)
    inter_community_edges_configs = 0
    mue_seq = [Decimal(mue) for mue in range(phi)]
    #for mue in mue_seq:
    if True:
        inter_community_edges_configs += comb(phi, mue) * (1-h_11) ** (mue) * (1-h_22) ** (phi - mue)



    stub_kind_1 = Decimal( stub_kind_1 )
    stub_kind_2 = Decimal( stub_kind_2 )

    pairing_prob =\
   total_pairing_num( stub_kind_1 - phi)\
 * total_pairing_num( stub_kind_2 - phi)\
 * comb( stub_kind_1, phi)\
 * comb( stub_kind_2, phi)\
 * factorial(phi)\
 * (2*h_11) ** ((stub_kind_1 - phi)/2)\
 * (2*h_22) ** ((stub_kind_2 - phi)/2)\
 * inter_community_edges_configs


#     h = Decimal(h)
#     pairing_prob = pairing_num * ( (1-h) / h ) ** phi
    return pairing_prob

def pairing_cond_num_asym_asym_mue_conditioned_norm_fact( stub_kind_1, stub_kind_2, h_11, h_22 ): #phi: number of heterogeneous pairs
    s_minority = np.min([stub_kind_1, stub_kind_2])

    min_phi = stub_kind_1 % 2 # if stub groups are odd, we need at least 1 pairing edge, otherwise 0 is the minimum

    phi_seq = range(min_phi, s_minority, 2)

    norm_factor = 0

    for phi in phi_seq:
        pairing_prob = pairing_cond_num_asym_asym_h( stub_kind_1, stub_kind_2, phi, h_11, h_22 )
#         print(pairing_prob)
        norm_factor += pairing_prob


    return norm_factor



def directed_asymmetric( s_1, s_2, m_12, m_21, h_11, h_22 ):

# why is this function so wrong?
# the correct version is currently titled: "pairing_cond_num_asym_asym_h_mue_conditioned"
#     print('m_21=',m_21)
#     print('m_12=',m_12)
    h_11 = Decimal(h_11)
    h_22 = Decimal(h_22)

    m_11 = (s_1 - m_21 - m_12)
    m_22 = (s_2 - m_21 - m_12)

    if ((m_22 % 2 != 0) or (m_11 % 2 != 0)):
        print('s and m values do not have suitable parity.')
        return 0
    elif (m_22 < 0) or (m_11 < 0):
        print('interlinks are higher than total stubs.')
        return 0

    else:
#         print('m_21=',m_21)
        m_11 = Decimal(m_11/2)
        m_22 = Decimal(m_22/2)

        h_12 = 1 - h_11
        h_21 = 1 - h_22
        # print(float(h_11), float(h_12), float(h_21), float(h_22))
        numerator = (h_11 ** m_11) * (h_12 ** m_12) * (h_21 ** m_21) * (h_22 ** m_22)
        denominator = factorial(m_11) * factorial(m_12)* factorial(m_21)* factorial(m_22)
#         print('m_21=',m_21)
        # print (m_12, m_21, ' num=', numerator,' den=', denominator)
#         print ((h_11 , m_11) , (h_12 , m_12) , (h_21 , m_21) , (h_22 , m_22))
        # print (('m_11\t' ) , ('m_12\t' ) , ('m_21\t' ) , ('m_22\t' ))
        # print ((m_11 ),'\t' , (m_12 ),'\t' , (m_21 ),'\t' , (m_22 ))
        # print('\n')
        return numerator/denominator

def directed_asymmetric_norm_fact( s_1, s_2, h_11, h_22 ): #phi: number of heterogeneous pairs
    if ((s_1 + s_2) % 2 == 0):
        min_m = s_1 % 2
        s_minority = np.min([s_1, s_2])
        m_seq = range(min_m, s_minority+1, 2)

        norm_factor = 0

        for m in m_seq:
            for m_12 in range(m+1):
                m_21 = m - m_12
                pairing_prob = directed_asymmetric( s_1, s_2, m_12, m_21, h_11, h_22 )
#                 print('m21','m12', m_21,m_12)
                norm_factor += pairing_prob
#                 print( s_1, s_2, m_12, m_21, h_11, h_22, pairing_prob )

        print(norm_factor)
        return norm_factor


    else:
        print('s values are not valid.')
        return 0

def directed_asymmetric_growth( s_1_o, s_2_o, m_12, m_21, h_11, h_22, k, k2 = None ): #phi: number of heterogeneous pairs
    if ((s_1_o % k != 0) or (s_2_o % k != 0)):
        print('invalid outgoing stubs')
        return 0
    elif (m_12 > s_1_o) or (m_21 > s_2_o):
        print('invalid outgoing stubs')
        return 0

    else:
        h_11 = Decimal(h_11)
        h_22 = Decimal(h_22)
        n_1 = Decimal(s_1_o / k)
        if k2:
            n_2 = Decimal(s_2_o / k2)
        else:
            n_2 = Decimal(s_2_o / k)
        h_21 = Decimal(1 - h_22)
        h_12 = Decimal(1 - h_11)
        rho_11 = Decimal((n_1 * h_11) / ((n_1 * h_11) + n_2 * h_12))
        rho_22 = Decimal((n_2 * h_22) / ((n_2 * h_22) + n_2 * h_21))
        rho_12 = Decimal(1 - rho_11)
        rho_21 = Decimal(1 - rho_22)

        m_11 = s_1_o - m_12
        m_22 = s_2_o - m_21

        probability = comb(s_1_o, m_12)\
        *comb(s_2_o, m_21)\
        *rho_11 ** m_11\
        *rho_12 ** m_12\
        *rho_21 ** m_21\
        *rho_22 ** m_22

        return probability

def directed_asymmetric_growth_norm_fact( s_1_o, s_2_o, h_11, h_22, k, k2 = None ): #phi: number of heterogeneous pairs
    norm_fact = 0
    for m_12 in range(s_1_o + 1):
        for m_21 in range(s_2_o + 1):
            prob = directed_asymmetric_growth( s_1_o, s_2_o, m_12, m_21, h_11, h_22, k , k2)
#             print(prob)
            norm_fact += prob
    return norm_fact


def log_pairing_num_asym_max( stub_kind_1, stub_kind_2, h ): #phi: number of heterogeneous pairs
    s_minority = np.min([stub_kind_1, stub_kind_2])
    phi_seq = range(0, s_minority, 2)


    log_pairing_prob = np.zeros_like(phi_seq, dtype = float)

    for i, phi in enumerate(phi_seq):
#         print(phi)
        pairing_num = log_total_pairing_num( stub_kind_1 - phi) \
     + log_total_pairing_num( stub_kind_2 - phi)\
     + log_comb( stub_kind_1, phi)\
     + log_comb( stub_kind_2, phi)\
     + log_factorial(phi)

        pairing_prob = pairing_num + np.log( (1-h) / h ) * phi
#         print(pairing_prob)
#         log_norm_factor += pairing_prob
        log_pairing_prob[ i ] = pairing_prob
#     print( log_pairing_prob )
    max_i = log_pairing_prob.argmax()

    return float(phi_seq[max_i] / s_minority)

##Hypergraphs:
def triparing_num( N_0, N_1, t_0, t_3, h ):
    N_0 = Decimal(N_0)
    N_1 = Decimal(N_1)
    t_0 = Decimal(t_0)
    t_3 = Decimal(t_3)
    # if h != 1: #don't know why!
    if True:
        h = Decimal(h)
    mixed_tris = Decimal( N_0 + N_1 - 3 * t_0 - 3 * t_3 )


    r_0 = N_0 - 3 * t_0
    r_1 = N_1 - 3 * t_3

    t_1 = (2*r_0 - r_1) / 3
    t_2 = (2*r_1 - r_0) / 3
    # print(mixed_tris == r_0 + r_1)


    # print(t_1, t_2)
    if (t_1 < 0) or (t_2 < 0):
        return 0
    if ((t_1 % 1) != 0) or ((t_2 % 1) != 0):
        print('remainder issue!')

#     t_2 = Decimal(t_2)
#     t_1 = Decimal(t_1)

#     r_0 = Decimal(r_0)
#     r_1 = Decimal(r_1)

    # print(t_0, t_1, t_2, t_3, N_0, N_1, h)
    # print('t: ',t_0, t_1, 'r: ', r_0, r_1)

    if(mixed_tris != r_0 + r_1):
        print('mixed_tris counted wrong')
    return\
      comb(N_0, 3 * t_0)\
    * comb(N_1, 3 * t_3)\
    * factorial(3 * t_0)\
    * factorial(3 * t_3)\
    / factorial(t_0)\
    / factorial(t_3)\
    / (6)**(t_0)\
    / (6)**(t_3)\
    * h ** (t_3 + t_0)\
    * (1-h) ** mixed_tris\
    * comb(r_1, t_1)\
    * comb(r_0, t_2)\
    * factorial( 2 * t_1 )\
    * 2 ** t_1\
    * factorial( 2 * t_2 )\
    * 2 ** t_2

def triparing_num_norm(N_0, N_1,h):
    max_t_3 = N_1 // 3
    max_t_0 = N_0 // 3
    norm_fact = 0
    for t_0 in range(0, max_t_0 + 1):

        left_over_0 = N_0 - (t_0 * 3)
        min_t_3 = np.max([0, (N_1 - 2*left_over_0)// 3])
        # max_t_3 = (N_1 - left_over_0)// 3
        # for t_3 in range(min_t_3, max_t_3 + 1):
        for t_3 in range(0, max_t_3 + 1):
            # print(t_0, t_3)
            non_norm_prob = triparing_num( N_0, N_1, t_0, t_3, h )
            norm_fact += non_norm_prob
    return norm_fact

if __name__ == "__main__":
    #print( log_pairing_cond_num_asym(40,40,10,0.5) )
    #stub_kind_1, stub_kind_2, phi, h = 100000, 100000, 50, 0.5
    #print( pairing_cond_num_asym( float(stub_kind_1), stub_kind_2, phi, h) )
    #total_pairing_num(100)
    #print('hi')
    print( pairing_cond_num_asym_norm_fact(20, 20, 5) )
