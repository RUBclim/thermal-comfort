"""Adjusted code by Fredrik Lindberg, orginally from Dag Wastberg, Rambolls"""
from collections import namedtuple

import numpy as np
from numba import njit

PET_person = namedtuple('PET_person', 'mbody age height activity sex clo')


def calculate_PET_index(Ta, Pa, Tmrt, va, pet):
    mbody = pet.mbody
    age = pet.age
    height = pet.height
    activity = pet.activity
    sex = pet.sex
    clo = pet.clo
    pet_index = np.zeros_like(Tmrt)
    for y in range(pet_index.shape[0]):
        for x in range(pet_index.shape[1]):
            pet_index[y, y] = _PET(
                Ta[x], Pa[x], Tmrt[x][y], va[x]
                [y], mbody, age, height, activity, clo, sex,
            )


def calculate_PET_grid(Ta, RH, Tmrt, va, pet):
    mbody = pet.mbody
    age = pet.age
    height = pet.height
    activity = pet.activity
    sex = pet.sex
    clo = pet.clo
    pet_index = np.zeros_like(Tmrt)
    total = 100. / (int(pet_index.shape[0] * pet_index.shape[1]))
    index = 0
    # print(Tmrt.shape)
    # print(va.shape)
    for y in range(pet_index.shape[0]):
        for x in range(pet_index.shape[1]):
            if va[y, x] > 0:
                index = index + 1
                pet_index[y, x] = _PET(
                    Ta, RH, Tmrt[y, x], va[y, x],
                    mbody, age, height, activity, clo, sex,
                )
            else:
                pet_index[y, x] = -9999

    return pet_index


def calculate_PET_index_vec(Ta, Pa, Tmrt, va, pet):
    mbody = pet.mbody
    age = pet.age
    height = pet.height
    activity = pet.activity
    sex = pet.sex
    clo = pet.clo

    pet_index = _PET(Ta, Pa, Tmrt, va, mbody, age, height, activity, clo, sex)


@njit
def _PET(ta, RH, tmrt, v, mbody, age, ht, work, icl, sex):
    """
    Args:
        ta: air temperature
        RH: relative humidity
        tmrt: Mean Radiant temperature
        v: wind at pedestrian heigh
        mbody: body masss (kg)
        age: person's age (years)
        ht: height (meters)
        work: activity level (W)
        icl: clothing amount (0-5)
        sex: 1=male 2=female
    Returns:
    """

    # humidity conversion
    vps = 6.107 * (10. ** (7.5 * ta / (238. + ta)))
    vpa = RH * vps / 100  # water vapour presure, kPa

    po = 1013.25  # Pressure
    p = 1013.25  # Pressure
    rob = 1.06
    cb = 3.64 * 1000
    food = 0
    emsk = 0.99
    emcl = 0.95
    evap = 2.42e6
    sigma = 5.67e-8
    cair = 1.01 * 1000

    eta = 0  # No idea what eta is

    c_1 = 0.
    c_2 = 0.
    c_3 = 0.
    c_4 = 0.
    c_5 = 0.
    c_6 = 0.
    c_7 = 0.
    c_8 = 0.
    c_9 = 0.
    c_10 = 0.
    c_11 = 0.

    # INBODY
    metbf = 3.19 * mbody ** (3 / 4) * (
        1 + 0.004 * (30 - age) +
        0.018 * ((ht * 100 / (mbody ** (1 / 3))) - 42.1)
    )
    metbm = 3.45 * mbody ** (3 / 4) * (
        1 + 0.004 * (30 - age) +
        0.010 * ((ht * 100 / (mbody ** (1 / 3))) - 43.4)
    )
    if sex == 1:
        met = metbm + work
    else:
        met = metbf + work

    h = met * (1 - eta)
    rtv = 1.44e-6 * met

    # sensible respiration energy
    tex = 0.47 * ta + 21.0
    eres = cair * (ta - tex) * rtv

    # latent respiration energy
    vpex = 6.11 * 10 ** (7.45 * tex / (235 + tex))
    erel = 0.623 * evap / p * (vpa - vpex) * rtv
    # sum of the results
    ere = eres + erel

    # calcul constants
    feff = 0.725
    adu = 0.203 * mbody ** 0.425 * ht ** 0.725
    facl = (-2.36 + 173.51 * icl - 100.76 * icl * icl + 19.28 * (icl ** 3)) / 100
    if facl > 1:
        facl = 1
    rcl = (icl / 6.45) / facl
    y = 1

    # should these be else if statements?
    if icl < 2:
        y = (ht-0.2) / ht
    if icl <= 0.6:
        y = 0.5
    if icl <= 0.3:
        y = 0.1

    fcl = 1 + 0.15 * icl
    r2 = adu * (fcl - 1. + facl) / (2 * 3.14 * ht * y)
    r1 = facl * adu / (2 * 3.14 * ht * y)
    di = r2 - r1
    acl = adu * facl + adu * (fcl - 1)

    tcore = [0] * 8

    wetsk = 0
    hc = 2.67 + 6.5 * v ** 0.67
    hc = hc * (p / po) ** 0.55
    c_1 = h + ere
    he = 0.633 * hc / (p * cair)
    fec = 1 / (1 + 0.92 * hc * rcl)
    htcl = 6.28 * ht * y * di / (rcl * np.log(r2 / r1) * acl)
    aeff = adu * feff
    c_2 = adu * rob * cb
    c_5 = 0.0208 * c_2
    c_6 = 0.76075 * c_2
    rdsk = 0.79 * 10 ** 7
    rdcl = 0

    count2 = 0
    j = 1

    while count2 == 0 and j < 7:
        tsk = 34
        count1 = 0
        tcl = (ta + tmrt + tsk) / 3
        count3 = 1
        enbal2 = 0

        while count1 <= 3:
            enbal = 0
            while (enbal*enbal2) >= 0 and count3 < 200:
                enbal2 = enbal
                # 20
                rclo2 = emcl * sigma * ((tcl + 273.2) ** 4 - (tmrt + 273.2) ** 4) * feff
                tsk = 1 / htcl * (hc * (tcl - ta) + rclo2) + tcl

                # radiation balance
                rbare = aeff * (1 - facl) * emsk * sigma * \
                    ((tmrt + 273.2) ** 4 - (tsk + 273.2) ** 4)
                rclo = feff * acl * emcl * sigma * \
                    ((tmrt + 273.2) ** 4 - (tcl + 273.2) ** 4)
                rsum = rbare + rclo

                # convection
                cbare = hc * (ta - tsk) * adu * (1 - facl)
                cclo = hc * (ta - tcl) * acl
                csum = cbare + cclo

                # core temperature
                c_3 = 18 - 0.5 * tsk
                c_4 = 5.28 * adu * c_3
                c_7 = c_4 - c_6 - tsk * c_5
                c_8 = -c_1 * c_3 - tsk * c_4 + tsk * c_6
                c_9 = c_7 * c_7 - 4. * c_5 * c_8
                c_10 = 5.28 * adu - c_6 - c_5 * tsk
                c_11 = c_10 * c_10 - 4 * c_5 * (c_6 * tsk - c_1 - 5.28 * adu * tsk)
                # tsk[tsk==36]=36.01
                # print(tsk.shape)
                if tsk == 36:
                    tsk = 36.01

                tcore[7] = c_1 / (5.28 * adu + c_2 * 6.3 / 3600) + tsk
                tcore[3] = c_1 / (
                    5.28 * adu + (c_2 * 6.3 / 3600) /
                    (1 + 0.5 * (34 - tsk))
                ) + tsk
                if c_11 >= 0:
                    tcore[6] = (-c_10-c_11 ** 0.5) / (2 * c_5)
                if c_11 >= 0:
                    tcore[1] = (-c_10+c_11 ** 0.5) / (2 * c_5)
                if c_9 >= 0:
                    tcore[2] = (-c_7+abs(c_9) ** 0.5) / (2 * c_5)
                if c_9 >= 0:
                    tcore[5] = (-c_7-abs(c_9) ** 0.5) / (2 * c_5)
                tcore[4] = c_1 / (5.28 * adu + c_2 * 1 / 40) + tsk

                # transpiration
                tbody = 0.1 * tsk + 0.9 * tcore[j]
                sw = 304.94 * (tbody - 36.6) * adu / 3600000
                vpts = 6.11 * 10 ** (7.45 * tsk / (235. + tsk))
                if tbody <= 36.6:
                    sw = 0
                if sex == 2:
                    sw = 0.7 * sw
                eswphy = -sw * evap

                eswpot = he * (vpa - vpts) * adu * evap * fec
                wetsk = eswphy / eswpot
                if wetsk > 1:
                    wetsk = 1
                eswdif = eswphy - eswpot
                if eswdif <= 0:
                    esw = eswpot
                else:
                    esw = eswphy
                if esw > 0:
                    esw = 0

                # diffusion
                ed = evap / (rdsk + rdcl) * adu * (1 - wetsk) * (vpa - vpts)

                # MAX VB
                vb1 = 34 - tsk
                vb2 = tcore[j] - 36.6
                if vb2 < 0:
                    vb2 = 0
                if vb1 < 0:
                    vb1 = 0
                vb = (6.3 + 75 * vb2) / (1 + 0.5 * vb1)

                # energy balance
                enbal = h + ed + ere + esw + csum + rsum + food

                # clothing's temperature
                if count1 == 0:
                    xx = 1
                if count1 == 1:
                    xx = 0.1
                if count1 == 2:
                    xx = 0.01
                if count1 == 3:
                    xx = 0.001
                if enbal > 0:
                    tcl = tcl + xx
                else:
                    tcl = tcl - xx

                count3 = count3 + 1
            count1 = count1 + 1
            enbal2 = 0

        if j == 2 or j == 5:
            if c_9 >= 0:
                if tcore[j] >= 36.6 and tsk <= 34.050:
                    if (j != 4 and vb >= 91) or (j == 4 and vb < 89):
                        pass
                    else:
                        if vb > 90:
                            vb = 90
                        count2 = 1

        if j == 6 or j == 1:
            if c_11 > 0:
                if tcore[j] >= 36.6 and tsk > 33.850:
                    if (j != 4 and vb >= 91) or (j == 4 and vb < 89):
                        pass
                    else:
                        if vb > 90:
                            vb = 90
                        count2 = 1

        if j == 3:
            if tcore[j] < 36.6 and tsk <= 34.000:
                if (j != 4 and vb >= 91) or (j == 4 and vb < 89):
                    pass
                else:
                    if vb > 90:
                        vb = 90
                    count2 = 1

        if j == 7:
            if tcore[j] < 36.6 and tsk > 34.000:
                if (j != 4 and vb >= 91) or (j == 4 and vb < 89):
                    pass
                else:
                    if vb > 90:
                        vb = 90
                    count2 = 1

        if j == 4:
            if (j != 4 and vb >= 91) or (j == 4 and vb < 89):
                pass
            else:
                if vb > 90:
                    vb = 90
                count2 = 1

        j = j + 1

    # PET_cal
    tx = ta
    enbal2 = 0
    count1 = 0
    enbal = 0

    hc = 2.67 + 6.5 * 0.1 ** 0.67
    hc = hc * (p / po) ** 0.55

    while count1 <= 3:
        while (enbal * enbal2) >= 0:
            enbal2 = enbal

            # radiation balance
            rbare = aeff * (1 - facl) * emsk * sigma * \
                ((tx + 273.2) ** 4 - (tsk + 273.2) ** 4)
            rclo = feff * acl * emcl * sigma * ((tx + 273.2) ** 4 - (tcl + 273.2) ** 4)
            rsum = rbare + rclo

            # convection
            cbare = hc * (tx - tsk) * adu * (1 - facl)
            cclo = hc * (tx - tcl) * acl
            csum = cbare + cclo

            # diffusion
            ed = evap / (rdsk + rdcl) * adu * (1 - wetsk) * (12 - vpts)

            # respiration
            tex = 0.47 * tx + 21
            eres = cair * (tx - tex) * rtv
            vpex = 6.11 * 10 ** (7.45 * tex / (235 + tex))
            erel = 0.623 * evap / p * (12 - vpex) * rtv
            ere = eres + erel

            # energy balance
            enbal = h + ed + ere + esw + csum + rsum

            # iteration concerning Tx
            if count1 == 0:
                xx = 1
            if count1 == 1:
                xx = 0.1
            if count1 == 2:
                xx = 0.01
            if count1 == 3:
                xx = 0.001
            if enbal > 0:
                tx = tx - xx
            if enbal < 0:
                tx = tx + xx
        count1 = count1 + 1
        enbal2 = 0

    return tx


def utci_polynomial(D_Tmrt, Ta, va, Pa):

    # calculate 6th order polynomial as approximation
    UTCI_approx = Ta + \
        (6.07562052E-01) + \
        (-2.27712343E-02) * Ta + \
        (8.06470249E-04) * Ta * Ta + \
        (-1.54271372E-04) * Ta * Ta * Ta + \
        (-3.24651735E-06) * Ta * Ta * Ta * Ta + \
        (7.32602852E-08) * Ta * Ta * Ta * Ta * Ta + \
        (1.35959073E-09) * Ta * Ta * Ta * Ta * Ta * Ta + \
        (-2.25836520E+00) * va + \
        (8.80326035E-02) * Ta * va + \
        (2.16844454E-03) * Ta * Ta * va + \
        (-1.53347087E-05) * Ta * Ta * Ta * va + \
        (-5.72983704E-07) * Ta * Ta * Ta * Ta * va + \
        (-2.55090145E-09) * Ta * Ta * Ta * Ta * Ta * va + \
        (-7.51269505E-01) * va * va + \
        (-4.08350271E-03) * Ta * va * va + \
        (-5.21670675E-05) * Ta * Ta * va * va + \
        (1.94544667E-06) * Ta * Ta * Ta * va * va + \
        (1.14099531E-08) * Ta * Ta * Ta * Ta * va * va + \
        (1.58137256E-01) * va * va * va + \
        (-6.57263143E-05) * Ta * va * va * va + \
        (2.22697524E-07) * Ta * Ta * va * va * va + \
        (-4.16117031E-08) * Ta * Ta * Ta * va * va * va + \
        (-1.27762753E-02) * va * va * va * va + \
        (9.66891875E-06) * Ta * va * va * va * va + \
        (2.52785852E-09) * Ta * Ta * va * va * va * va + \
        (4.56306672E-04) * va * va * va * va * va + \
        (-1.74202546E-07) * Ta * va * va * va * va * va + \
        (-5.91491269E-06) * va * va * va * va * va * va + \
        (3.98374029E-01) * D_Tmrt + \
        (1.83945314E-04) * Ta * D_Tmrt + \
        (-1.73754510E-04) * Ta * Ta * D_Tmrt + \
        (-7.60781159E-07) * Ta * Ta * Ta * D_Tmrt + \
        (3.77830287E-08) * Ta * Ta * Ta * Ta * D_Tmrt + \
        (5.43079673E-10) * Ta * Ta * Ta * Ta * Ta * D_Tmrt + \
        (-2.00518269E-02) * va * D_Tmrt + \
        (8.92859837E-04) * Ta * va * D_Tmrt + \
        (3.45433048E-06) * Ta * Ta * va * D_Tmrt + \
        (-3.77925774E-07) * Ta * Ta * Ta * va * D_Tmrt + \
        (-1.69699377E-09) * Ta * Ta * Ta * Ta * va * D_Tmrt + \
        (1.69992415E-04) * va * va * D_Tmrt + \
        (-4.99204314E-05) * Ta * va * va * D_Tmrt + \
        (2.47417178E-07) * Ta * Ta * va * va * D_Tmrt + \
        (1.07596466E-08) * Ta * Ta * Ta * va * va * D_Tmrt + \
        (8.49242932E-05) * va * va * va * D_Tmrt + \
        (1.35191328E-06) * Ta * va * va * va * D_Tmrt + \
        (-6.21531254E-09) * Ta * Ta * va * va * va * D_Tmrt + \
        (-4.99410301E-06) * va * va * va * va * D_Tmrt + \
        (-1.89489258E-08) * Ta * va * va * va * va * D_Tmrt + \
        (8.15300114E-08) * va * va * va * va * va * D_Tmrt + \
        (7.55043090E-04) * D_Tmrt * D_Tmrt + \
        (-5.65095215E-05) * Ta * D_Tmrt * D_Tmrt + \
        (-4.52166564E-07) * Ta * Ta * D_Tmrt * D_Tmrt + \
        (2.46688878E-08) * Ta * Ta * Ta * D_Tmrt * D_Tmrt + \
        (2.42674348E-10) * Ta * Ta * Ta * Ta * D_Tmrt * D_Tmrt + \
        (1.54547250E-04) * va * D_Tmrt * D_Tmrt + \
        (5.24110970E-06) * Ta * va * D_Tmrt * D_Tmrt + \
        (-8.75874982E-08) * Ta * Ta * va * D_Tmrt * D_Tmrt + \
        (-1.50743064E-09) * Ta * Ta * Ta * va * D_Tmrt * D_Tmrt + \
        (-1.56236307E-05) * va * va * D_Tmrt * D_Tmrt + \
        (-1.33895614E-07) * Ta * va * va * D_Tmrt * D_Tmrt + \
        (2.49709824E-09) * Ta * Ta * va * va * D_Tmrt * D_Tmrt + \
        (6.51711721E-07) * va * va * va * D_Tmrt * D_Tmrt + \
        (1.94960053E-09) * Ta * va * va * va * D_Tmrt * D_Tmrt + \
        (-1.00361113E-08) * va * va * va * va * D_Tmrt * D_Tmrt + \
        (-1.21206673E-05) * D_Tmrt * D_Tmrt * D_Tmrt + \
        (-2.18203660E-07) * Ta * D_Tmrt * D_Tmrt * D_Tmrt + \
        (7.51269482E-09) * Ta * Ta * D_Tmrt * D_Tmrt * D_Tmrt + \
        (9.79063848E-11) * Ta * Ta * Ta * D_Tmrt * D_Tmrt * D_Tmrt + \
        (1.25006734E-06) * va * D_Tmrt * D_Tmrt * D_Tmrt + \
        (-1.81584736E-09) * Ta * va * D_Tmrt * D_Tmrt * D_Tmrt + \
        (-3.52197671E-10) * Ta * Ta * va * D_Tmrt * D_Tmrt * D_Tmrt + \
        (-3.36514630E-08) * va * va * D_Tmrt * D_Tmrt * D_Tmrt + \
        (1.35908359E-10) * Ta * va * va * D_Tmrt * D_Tmrt * D_Tmrt + \
        (4.17032620E-10) * va * va * va * D_Tmrt * D_Tmrt * D_Tmrt + \
        (-1.30369025E-09) * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt + \
        (4.13908461E-10) * Ta * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt + \
        (9.22652254E-12) * Ta * Ta * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt + \
        (-5.08220384E-09) * va * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt + \
        (-2.24730961E-11) * Ta * va * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt + \
        (1.17139133E-10) * va * va * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt + \
        (6.62154879E-10) * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt + \
        (4.03863260E-13) * Ta * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt + \
        (1.95087203E-12) * va * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt + \
        (-4.73602469E-12) * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt + \
        (5.12733497E+00) * Pa + \
        (-3.12788561E-01) * Ta * Pa + \
        (-1.96701861E-02) * Ta * Ta * Pa + \
        (9.99690870E-04) * Ta * Ta * Ta * Pa + \
        (9.51738512E-06) * Ta * Ta * Ta * Ta * Pa + \
        (-4.66426341E-07) * Ta * Ta * Ta * Ta * Ta * Pa + \
        (5.48050612E-01) * va * Pa + \
        (-3.30552823E-03) * Ta * va * Pa + \
        (-1.64119440E-03) * Ta * Ta * va * Pa + \
        (-5.16670694E-06) * Ta * Ta * Ta * va * Pa + \
        (9.52692432E-07) * Ta * Ta * Ta * Ta * va * Pa + \
        (-4.29223622E-02) * va * va * Pa + \
        (5.00845667E-03) * Ta * va * va * Pa + \
        (1.00601257E-06) * Ta * Ta * va * va * Pa + \
        (-1.81748644E-06) * Ta * Ta * Ta * va * va * Pa + \
        (-1.25813502E-03) * va * va * va * Pa + \
        (-1.79330391E-04) * Ta * va * va * va * Pa + \
        (2.34994441E-06) * Ta * Ta * va * va * va * Pa + \
        (1.29735808E-04) * va * va * va * va * Pa + \
        (1.29064870E-06) * Ta * va * va * va * va * Pa + \
        (-2.28558686E-06) * va * va * va * va * va * Pa + \
        (-3.69476348E-02) * D_Tmrt * Pa + \
        (1.62325322E-03) * Ta * D_Tmrt * Pa + \
        (-3.14279680E-05) * Ta * Ta * D_Tmrt * Pa + \
        (2.59835559E-06) * Ta * Ta * Ta * D_Tmrt * Pa + \
        (-4.77136523E-08) * Ta * Ta * Ta * Ta * D_Tmrt * Pa + \
        (8.64203390E-03) * va * D_Tmrt * Pa + \
        (-6.87405181E-04) * Ta * va * D_Tmrt * Pa + \
        (-9.13863872E-06) * Ta * Ta * va * D_Tmrt * Pa + \
        (5.15916806E-07) * Ta * Ta * Ta * va * D_Tmrt * Pa + \
        (-3.59217476E-05) * va * va * D_Tmrt * Pa + \
        (3.28696511E-05) * Ta * va * va * D_Tmrt * Pa + \
        (-7.10542454E-07) * Ta * Ta * va * va * D_Tmrt * Pa + \
        (-1.24382300E-05) * va * va * va * D_Tmrt * Pa + \
        (-7.38584400E-09) * Ta * va * va * va * D_Tmrt * Pa + \
        (2.20609296E-07) * va * va * va * va * D_Tmrt * Pa + \
        (-7.32469180E-04) * D_Tmrt * D_Tmrt * Pa + \
        (-1.87381964E-05) * Ta * D_Tmrt * D_Tmrt * Pa + \
        (4.80925239E-06) * Ta * Ta * D_Tmrt * D_Tmrt * Pa + \
        (-8.75492040E-08) * Ta * Ta * Ta * D_Tmrt * D_Tmrt * Pa + \
        (2.77862930E-05) * va * D_Tmrt * D_Tmrt * Pa + \
        (-5.06004592E-06) * Ta * va * D_Tmrt * D_Tmrt * Pa + \
        (1.14325367E-07) * Ta * Ta * va * D_Tmrt * D_Tmrt * Pa + \
        (2.53016723E-06) * va * va * D_Tmrt * D_Tmrt * Pa + \
        (-1.72857035E-08) * Ta * va * va * D_Tmrt * D_Tmrt * Pa + \
        (-3.95079398E-08) * va * va * va * D_Tmrt * D_Tmrt * Pa + \
        (-3.59413173E-07) * D_Tmrt * D_Tmrt * D_Tmrt * Pa + \
        (7.04388046E-07) * Ta * D_Tmrt * D_Tmrt * D_Tmrt * Pa + \
        (-1.89309167E-08) * Ta * Ta * D_Tmrt * D_Tmrt * D_Tmrt * Pa + \
        (-4.79768731E-07) * va * D_Tmrt * D_Tmrt * D_Tmrt * Pa + \
        (7.96079978E-09) * Ta * va * D_Tmrt * D_Tmrt * D_Tmrt * Pa + \
        (1.62897058E-09) * va * va * D_Tmrt * D_Tmrt * D_Tmrt * Pa + \
        (3.94367674E-08) * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * Pa + \
        (-1.18566247E-09) * Ta * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * Pa + \
        (3.34678041E-10) * va * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * Pa + \
        (-1.15606447E-10) * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * Pa + \
        (-2.80626406E+00) * Pa * Pa + \
        (5.48712484E-01) * Ta * Pa * Pa + \
        (-3.99428410E-03) * Ta * Ta * Pa * Pa + \
        (-9.54009191E-04) * Ta * Ta * Ta * Pa * Pa + \
        (1.93090978E-05) * Ta * Ta * Ta * Ta * Pa * Pa + \
        (-3.08806365E-01) * va * Pa * Pa + \
        (1.16952364E-02) * Ta * va * Pa * Pa + \
        (4.95271903E-04) * Ta * Ta * va * Pa * Pa + \
        (-1.90710882E-05) * Ta * Ta * Ta * va * Pa * Pa + \
        (2.10787756E-03) * va * va * Pa * Pa + \
        (-6.98445738E-04) * Ta * va * va * Pa * Pa + \
        (2.30109073E-05) * Ta * Ta * va * va * Pa * Pa + \
        (4.17856590E-04) * va * va * va * Pa * Pa + \
        (-1.27043871E-05) * Ta * va * va * va * Pa * Pa + \
        (-3.04620472E-06) * va * va * va * va * Pa * Pa + \
        (5.14507424E-02) * D_Tmrt * Pa * Pa + \
        (-4.32510997E-03) * Ta * D_Tmrt * Pa * Pa + \
        (8.99281156E-05) * Ta * Ta * D_Tmrt * Pa * Pa + \
        (-7.14663943E-07) * Ta * Ta * Ta * D_Tmrt * Pa * Pa + \
        (-2.66016305E-04) * va * D_Tmrt * Pa * Pa + \
        (2.63789586E-04) * Ta * va * D_Tmrt * Pa * Pa + \
        (-7.01199003E-06) * Ta * Ta * va * D_Tmrt * Pa * Pa + \
        (-1.06823306E-04) * va * va * D_Tmrt * Pa * Pa + \
        (3.61341136E-06) * Ta * va * va * D_Tmrt * Pa * Pa + \
        (2.29748967E-07) * va * va * va * D_Tmrt * Pa * Pa + \
        (3.04788893E-04) * D_Tmrt * D_Tmrt * Pa * Pa + \
        (-6.42070836E-05) * Ta * D_Tmrt * D_Tmrt * Pa * Pa + \
        (1.16257971E-06) * Ta * Ta * D_Tmrt * D_Tmrt * Pa * Pa + \
        (7.68023384E-06) * va * D_Tmrt * D_Tmrt * Pa * Pa + \
        (-5.47446896E-07) * Ta * va * D_Tmrt * D_Tmrt * Pa * Pa + \
        (-3.59937910E-08) * va * va * D_Tmrt * D_Tmrt * Pa * Pa + \
        (-4.36497725E-06) * D_Tmrt * D_Tmrt * D_Tmrt * Pa * Pa + \
        (1.68737969E-07) * Ta * D_Tmrt * D_Tmrt * D_Tmrt * Pa * Pa + \
        (2.67489271E-08) * va * D_Tmrt * D_Tmrt * D_Tmrt * Pa * Pa + \
        (3.23926897E-09) * D_Tmrt * D_Tmrt * D_Tmrt * D_Tmrt * Pa * Pa + \
        (-3.53874123E-02) * Pa * Pa * Pa + \
        (-2.21201190E-01) * Ta * Pa * Pa * Pa + \
        (1.55126038E-02) * Ta * Ta * Pa * Pa * Pa + \
        (-2.63917279E-04) * Ta * Ta * Ta * Pa * Pa * Pa + \
        (4.53433455E-02) * va * Pa * Pa * Pa + \
        (-4.32943862E-03) * Ta * va * Pa * Pa * Pa + \
        (1.45389826E-04) * Ta * Ta * va * Pa * Pa * Pa + \
        (2.17508610E-04) * va * va * Pa * Pa * Pa + \
        (-6.66724702E-05) * Ta * va * va * Pa * Pa * Pa + \
        (3.33217140E-05) * va * va * va * Pa * Pa * Pa + \
        (-2.26921615E-03) * D_Tmrt * Pa * Pa * Pa + \
        (3.80261982E-04) * Ta * D_Tmrt * Pa * Pa * Pa + \
        (-5.45314314E-09) * Ta * Ta * D_Tmrt * Pa * Pa * Pa + \
        (-7.96355448E-04) * va * D_Tmrt * Pa * Pa * Pa + \
        (2.53458034E-05) * Ta * va * D_Tmrt * Pa * Pa * Pa + \
        (-6.31223658E-06) * va * va * D_Tmrt * Pa * Pa * Pa + \
        (3.02122035E-04) * D_Tmrt * D_Tmrt * Pa * Pa * Pa + \
        (-4.77403547E-06) * Ta * D_Tmrt * D_Tmrt * Pa * Pa * Pa + \
        (1.73825715E-06) * va * D_Tmrt * D_Tmrt * Pa * Pa * Pa + \
        (-4.09087898E-07) * D_Tmrt * D_Tmrt * D_Tmrt * Pa * Pa * Pa + \
        (6.14155345E-01) * Pa * Pa * Pa * Pa + \
        (-6.16755931E-02) * Ta * Pa * Pa * Pa * Pa + \
        (1.33374846E-03) * Ta * Ta * Pa * Pa * Pa * Pa + \
        (3.55375387E-03) * va * Pa * Pa * Pa * Pa + \
        (-5.13027851E-04) * Ta * va * Pa * Pa * Pa * Pa + \
        (1.02449757E-04) * va * va * Pa * Pa * Pa * Pa + \
        (-1.48526421E-03) * D_Tmrt * Pa * Pa * Pa * Pa + \
        (-4.11469183E-05) * Ta * D_Tmrt * Pa * Pa * Pa * Pa + \
        (-6.80434415E-06) * va * D_Tmrt * Pa * Pa * Pa * Pa + \
        (-9.77675906E-06) * D_Tmrt * D_Tmrt * Pa * Pa * Pa * Pa + \
        (8.82773108E-02) * Pa * Pa * Pa * Pa * Pa + \
        (-3.01859306E-03) * Ta * Pa * Pa * Pa * Pa * Pa + \
        (1.04452989E-03) * va * Pa * Pa * Pa * Pa * Pa + \
        (2.47090539E-04) * D_Tmrt * Pa * Pa * Pa * Pa * Pa + \
        (1.48348065E-03) * Pa * Pa * Pa * Pa * Pa * Pa

    return UTCI_approx


def utci_calculator(Ta, RH, Tmrt, va10m):
    # Program for calculating UTCI Temperature (UTCI)
    # released for public use after termination of COST Action 730

    # Translated from fortran by Fredrik Lindberg, Göteborg Urban Climate Group, Sweden
    # UTCI, Version a 0.002, October 2009
    # Copyright (C) 2009  Peter Broede

    if Ta <= -999 or RH <= -999 or va10m <= -999 or Tmrt <= -999:
        UTCI_approx = -999
    else:
        # saturation vapour pressure (es)
        g = np.array([
            -2.8365744E3, - 6.028076559E3, 1.954263612E1, - 2.737830188E-2,
            1.6261698E-5, 7.0229056E-10, - 1.8680009E-13, 2.7150305,
        ])

        tk = Ta + 273.15  # ! air temp in K
        es = g[7] * np.log(tk)
        for i in range(0, 7):
            es = es + g[i] * tk ** (i + 1 - 3.)

        es = np.exp(es) * 0.01

        ehPa = es * RH / 100.

        D_Tmrt = Tmrt - Ta
        Pa = ehPa / 10.0  # use vapour pressure in kPa
        va = va10m

        # Calculate 6th order polynomial as approximation
        UTCI_approx = utci_polynomial(D_Tmrt, Ta, va, Pa)

    return UTCI_approx


def utci_calculator_grid(Ta, RH, Tmrt, va10m):
    # Program for calculating UTCI Temperature (UTCI)
    # released for public use after termination of COST Action 730

    # Translated from fortran by Fredrik Lindberg, Göteborg Urban Climate Group, Sweden
    # UTCI, Version a 0.002, October 2009
    # Copyright (C) 2009  Peter Broede

    rows = Tmrt.shape[0]
    cols = Tmrt.shape[1]

    # If Tair or relative humidity is NaN, set all UTCI to -999
    if Ta <= -999 or RH <= -999:
        UTCI_approx = np.ones((rows, cols)) * -999
    # Else, create an empty raster for UTCI and start estimating UTCI
    else:
        UTCI_approx = np.zeros((rows, cols))

        # saturation vapour pressure (es)
        g = np.array([
            -2.8365744E3, - 6.028076559E3, 1.954263612E1, - 2.737830188E-2,
            1.6261698E-5, 7.0229056E-10, - 1.8680009E-13, 2.7150305,
        ])

        tk = Ta + 273.15  # ! air temp in K
        es = g[7] * np.log(tk)
        for i in range(0, 7):
            es = es + g[i] * tk ** (i + 1 - 3.)

        es = np.exp(es) * 0.01

        ehPa = es * RH / 100.
        Pa = ehPa / 10.0  # use vapour pressure in kPa

        # For progress counter
        index = 0
        total = 100. / (rows * cols)

        # Start UTCI approximation
        for iy in np.arange(rows):
            for ix in np.arange(cols):

                D_Tmrt = Tmrt[iy, ix] - Ta
                # va = va10m[iy, ix]

                # Set UTCI to NaN (-999) if Tmrt or wind speed is NaN
                if Tmrt[iy, ix] <= -999 or va10m[iy, ix] <= -999:
                    UTCI_approx[iy, ix] = -9999
                # Calculate 6th order polynomial as approximation if wind speed is above zero
                elif va10m[iy, ix] > 0:
                    UTCI_approx[iy, ix] = utci_polynomial(D_Tmrt, Ta, va10m[iy, ix], Pa)

                # Progress counter
                index = index + 1

    return UTCI_approx
