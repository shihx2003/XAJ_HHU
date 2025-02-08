# Decompiled with PyLingual (https://pylingual.io)
# Internal filename: XAJ_Water.py
# Bytecode version: 3.10.0rc2 (3439)
# Source timestamp: 1970-01-01 00:00:00 UTC (0)

def Yield(kc, ep, wu, wl, wd, c, b, wum, wlm, wm, w, imp, prec):
    """
    蒸散发 + 产流
    :param kc:
    :param ep:
    :param wu:
    :param wl:
    :param wd:
    :param c:
    :param b:
    :param wum:
    :param wlm:
    :param wm:
    :param w:
    :param EU:
    :param EL:
    :param ED:
    :param imp:
    :param prec:
    :param EK:
    :return:
    """
    wmm = (1 + b) * wm / (1 - imp)
    ek = ep * kc
    pe = prec - ek
    if abs(pe) < 0.001:
        pe = 0
    div = 5
    if pe > div:
        nd = int(pe / div) + 1
        ped = [pe / nd] * nd
    else:
        nd = 1
        ped = [pe / nd]
    rd = []
    if pe > 0:
        if abs(w - wm) < 0.001:
            a = wmm
        else:
            a = wmm * (1 - (1 - w / wm) ** (1 / (1 + b)))
        r = 0
        peds = 0
        for ind in range(nd):
            a = a + ped[ind]
            peds = peds + ped[ind]
            rds = r
            if a < wmm:
                r = peds - (wm - w) + wm * (1 - a / wmm) ** (1 + b)
            else:
                r = peds - (wm - w)
            rd.append(r - rds)
        eu = ek
        el = 0
        ed = 0
        if wu + pe - r > wum:
            if wu + pe - r - wum + wl > wlm:
                wu = wum
                wl = wlm
                wd = w + peds - r - wu - wl
            else:
                wl = wl + wu + pe - r - wum
                wu = wum
        else:
            wu = wu + pe - r
    else:
        r = 0
        rd.append(r)
        if wu + pe >= 0:
            eu = ek
            ed = 0
            el = 0
            wu = wu + pe
        else:
            eu = wu + ek + pe
            wu = 0
            if wl > c * wlm:
                el = (ek - eu) * wl / wlm
                ed = 0
                wl = wl - el
            elif wl > c * (ek - eu):
                el = c * (ek - eu)
                ed = 0
                wl = wl - el
            else:
                el = wl
                ed = c * (ek - eu) - el
                wl = 0
                wd = wd - ed
                if wd < 0:
                    wd = 0
    et = eu + el + ed
    w = wu + wl + wd
    r_yield = r
    return (et, eu, el, ed, w, wu, wl, wd, nd, ped, rd, r_yield)

def Divide3Source(IMP, SM, EX, KG, KI, PED, RD, ND, FR, S):
    """
    分水源计算

    :return:
    """
    PE = 0
    R = 0
    for iND1 in range(ND):
        PE = PE + PED[iND1]
        R = R + RD[iND1]
    if PE <= 0:
        RS = 0
        RG = S * KG * FR
        RI = S * KI * FR
        S = S * (1 - KG - KI)
    else:
        KID = (1 - (1 - (KG + KI)) ** (1 / ND)) / (KG + KI)
        KGD = KID * KG
        KID = KID * KI
        RS = 0
        RI = 0
        RG = 0
        SMM = (1 + EX) * SM
        RB = IMP * PE
        for iND2 in range(ND):
            TD = RD[iND2] - IMP * PED[iND2]
            X = FR
            FR = TD / PED[iND2]
            FR = max(1e-05, FR)
            S = X * S / FR
            RS_Adjust = 0
            if S > SM:
                RS_Adjust = (S - SM) * FR
                S = SM
            AU = SMM * (1 - (1 - S / SM) ** (1 / (1 + EX)))
            if AU + PED[iND2] < SMM:
                RSD = (PED[iND2] - SM + S + SM * (1 - (PED[iND2] + AU) / SMM) ** (1 + EX)) * FR
            else:
                RSD = (PED[iND2] + S - SM) * FR
            if RSD < 0:
                RSD = 0
            RS = RS + RSD + RS_Adjust
            S = S + PED[iND2] - RSD / FR
            RID = S * KID * FR
            RI = RI + RID
            RGD = S * KGD * FR
            RG = RG + RGD
            S = S + PED[iND2] - (RID + RGD) / FR
        RS = RS + RB
    return (RS, RI, RG, FR, S)

def hsRouting(CG, CI, RS, RI, RG, CP, QI, QG):
    """
    河网汇流
    :param ZoneParms:
    :param RS:
    :param RI:
    :param RG:
    :param CP:
    :param QS:
    :param QI:
    :param QG:
    :return:
    """
    QS = RS * CP
    QI = QI * CI + RI * (1 - CI) * CP
    QG = QG * CG + RG * (1 - CG) * CP
    return (QS, QI, QG)

def Chrouting(CS, LAG, QS, QI, QG, QS_CS, QI_CS, QG_CS):
    """
    河道汇流
    :param ZoneParms:
    :param QS:
    :param QI:
    :param QG:
    :param QS_CS:
    :param QI_CS:
    :return:
    """
    QS_CS = QS_CS * CS + QS * (1 - CS)
    QI_CS = QI_CS * CS + QI * (1 - CS)
    QG_CS = QG_CS * CS + QG * (1 - CS)
    return (QS_CS, QI_CS, QG_CS)

def Muskingum(dp, KE, XE, dt, O2, QXE):
    """
    马斯京根
    :param ZoneParms:
    :param dt:
    :param maxdp:
    :param O2:
    :param QXE:
    :return:
    """
    tempmp = dp + 1
    t = dt / 3600
    C0 = (0.5 * t - KE * XE) / (KE - KE * XE + 0.5 * t)
    C1 = (KE * XE + 0.5 * t) / (KE - KE * XE + 0.5 * t)
    C2 = (KE - KE * XE - 0.5 * t) / (KE - KE * XE + 0.5 * t)
    if tempmp == 1:
        QXE[tempmp - 1] = O2
        return (O2, QXE)
    for itm in range(2, tempmp + 1):
        itm = itm - 1
        I2 = O2
        I1 = QXE[itm - 1]
        O1 = QXE[itm]
        QXE[itm - 1] = O2
        O2 = C0 * I2 + C1 * I1 + C2 * O1
    QXE[tempmp - 1] = O2
    return (O2, QXE)