# Decompiled with PyLingual (https://pylingual.io)
# Internal filename: F:\Python_Project\HHU\Standard_XAJ\Test\XAJ_Driver.py
# Bytecode version: 3.10.0rc2 (3439)
# Source timestamp: 1970-01-01 00:00:00 UTC (0)

from pandas import read_csv, date_range, DataFrame
from numpy import array, zeros, tile, reshape
from f90nml import read
from datetime import datetime, timedelta
from XAJ_Water import *

from XAJ_Restart import *
from traceback import print_exc
work_dir = 'E:/Learning/XAJ model/XAJ_HHU/'
Prec_Dir = work_dir + 'Forcings/prec.csv'
Evapo_Dir = work_dir + 'Forcings/ep.csv'
Initials_Dir = work_dir + '/Parameters/parameters.csv'
Namelist_Dir = work_dir + '/Parameters/xaj.namelist'
Result_Dir = work_dir + '/Results/'
restart_dir = work_dir + '/Restarts/'
ini_p = 0.2
maxdp = 5

def Read_Inputs(Prec_Dir, Evapo_Dir):
    """
    Read the input of P, E
    Parameters:
    ============
    Input_fn: str, filename of input data including the location.

    Return:
    ============
    P, E, Q_obs: np.ndarray, precipitation, evaporation and observed flow discharge.
    """
    print('==============================Load the required Input Forcings=============================')
    Prec = read_csv(Prec_Dir)
    Prec = Prec.iloc[:, 1:]
    zones_columns = list(Prec.columns)
    Prec = array(Prec, dtype=float)
    Evapo = read_csv(Evapo_Dir)
    Evapo = Evapo.iloc[:, 1:]
    Evapo = array(Evapo, dtype=float)
    return (Prec, Evapo, zones_columns)

def read_parms(Initials_Dir):
    """
    Read runoff generation parameters and initial states.
    Parameters:
    ============
    f_parms: str, filename of parameters.
    f_initials: str, filename of initial states.

    Return:
    ============
    Values of parameters and initial states.
    """
    print('==================================Load the required parameters=================================')
    Initial_Parms = read_csv(Initials_Dir)
    print(Initial_Parms)
    AREA = array(Initial_Parms['Area'])
    DP = array(Initial_Parms['DP'])
    KC = array(Initial_Parms['KC'])
    B = array(Initial_Parms['B'])
    C = array(Initial_Parms['C'])
    IMP = array(Initial_Parms['IMP'])
    IMP[IMP > 0.5] = 0.5
    WM = array(Initial_Parms['WM'])
    WUM = array(Initial_Parms['WUM'])
    WLM = array(Initial_Parms['WLM'])
    WDM = WM - WUM - WLM
    SM = array(Initial_Parms['SM'])
    EX = array(Initial_Parms['EX'])
    KG = array(Initial_Parms['KG'])
    KI = array(Initial_Parms['KI'])
    CG = array(Initial_Parms['CG'])
    CI = array(Initial_Parms['CI'])
    CS = array(Initial_Parms['CS'])
    LAG = array(Initial_Parms['LAG'], dtype=int)
    KE = array(Initial_Parms['KE'])
    XE = array(Initial_Parms['XE'])
    return (AREA, DP, KC, B, C, IMP, WM, WUM, WLM, WDM, SM, EX, KG, KI, CG, CI, CS, LAG, KE, XE)

def XAJ_Namelist(Namelist_Dir):
    """
    Read Namelist.
    :param Namelist_Dir:
    :return:
    """
    print('===============================Load the required Namelist==============================')
    XAJ_namelist = read(Namelist_Dir)
    NZONE = XAJ_namelist['xaj_namelist']['nzone']
    print(f'nzone:{NZONE}')
    NSTEP = XAJ_namelist['xaj_namelist']['nstep']
    print(f'nstep:{NSTEP}')
    DT = XAJ_namelist['xaj_namelist']['dt']
    print(f'dt:{DT}')
    START_YEAR = XAJ_namelist['xaj_namelist']['start_year']
    print(f'year:{START_YEAR:04d}')
    START_MONTH = XAJ_namelist['xaj_namelist']['start_month']
    print(f'month:{START_MONTH:02d}')
    START_DAY = XAJ_namelist['xaj_namelist']['start_day']
    print(f'day:{START_DAY:02d}')
    START_HOUR = XAJ_namelist['xaj_namelist']['start_hour']
    print(f'hour:{START_HOUR:02d}')
    print(f'simulation start time:{START_YEAR:04d}4-{START_MONTH:02d}4-{START_DAY:02d}4 {START_HOUR:02d}4:00:00')
    try:
        restart_file = XAJ_namelist['xaj_namelist']['restart_file']
    except:
        restart_file = None
    DIR_RESTART_OUT = XAJ_namelist['xaj_namelist']['dir_restart_out']
    return (NZONE, NSTEP, DT, START_YEAR, START_MONTH, START_DAY, START_HOUR, restart_file, DIR_RESTART_OUT)




if __name__ == '__main__':
    try:
        NZONE, NSTEP, DT, START_YEAR, START_MONTH, START_DAY, START_HOUR, restart_file, DIR_RESTART_OUT = XAJ_Namelist(Namelist_Dir)
        p_area, p_dp, p_kc, p_b, p_c, p_imp, p_wm, p_wum, p_wlm, p_wdm, p_sm, p_ex, p_kg, p_ki, p_cg, p_ci, p_cs, p_lag, p_ke, p_xe = read_parms(Initials_Dir)
        run_datetime = datetime.now()
        Result_date = run_datetime.strftime('%Y%m%d%H%M%S')
        if DT != 86400:
            bb1 = p_ki + p_kg
            bb2 = p_ki / p_kg
            p_kg = (1 - (1 - bb1) ** (DT / 86400)) / (1 + bb2)
            p_ki = p_kg * bb2
            p_ci = p_ci ** (DT / 86400)
            p_cg = p_cg ** (DT / 86400)
        cp = 1000 * p_area / DT
        all_qg = []
        all_qig = []
        all_qsig = []
        if restart_file is None:
            start_date = str(START_YEAR).zfill(4) + '-' + str(START_MONTH).zfill(2) + '-' + str(START_DAY).zfill(2) + ' ' + str(START_HOUR).zfill(2) + ':' + str(0).zfill(2) + ':' + str(0).zfill(2)
            start_date = datetime.strptime(start_date, '%Y-%m-%d %H:%M:%S')
            end_date = datetime(START_YEAR, START_MONTH, START_DAY, START_HOUR) + timedelta(hours=float(NSTEP - 1))
            date_justification = date_range(start_date, end_date, freq='h')
            modeling_date_range = date_range(start_date, end_date, freq='h').astype(str)
            Prec, Evapo, zones_columns = Read_Inputs(Prec_Dir, Evapo_Dir)
            wu = ini_p * p_wum
            wl = ini_p * p_wlm
            wd = ini_p * p_wdm
            s = zeros(NZONE)
            fr = zeros(NZONE)
            qs = zeros(NZONE)
            qi = zeros(NZONE)
            qg = zeros(NZONE)
            qs_cs = zeros(NZONE)
            qi_cs = zeros(NZONE)
            qg_cs = zeros(NZONE)
            w = wu + wl + wd
            qs_lag = reshape(tile(zeros(NZONE), NSTEP + p_lag[0]), (NSTEP + p_lag[0], NZONE))
            qi_lag = reshape(tile(zeros(NZONE), NSTEP + p_lag[0]), (NSTEP + p_lag[0], NZONE))
            qg_lag = reshape(tile(zeros(NZONE), NSTEP + p_lag[0]), (NSTEP + p_lag[0], NZONE))
            qs_lag = array(qs_lag).T
            qi_lag = array(qi_lag).T
            qg_lag = array(qg_lag).T
            qxs_ini_nzone = reshape(tile(zeros(NZONE), maxdp), (maxdp, NZONE)).T
            qxi_ini_nzone = reshape(tile(zeros(NZONE), maxdp), (maxdp, NZONE)).T
            qxg_ini_nzone = reshape(tile(zeros(NZONE), maxdp), (maxdp, NZONE)).T
        else:
            Prec, Evapo, zones_columns, restart_date = get_restart_forcing(restart_file, Prec_Dir, Evapo_Dir)
            end_date = restart_date + timedelta(hours=float(NSTEP - 1))
            date_justification = date_range(restart_date, end_date, freq='h')
            modeling_date_range = date_range(restart_date, end_date, freq='h').astype(str)
            NSTEP = len(modeling_date_range)
            wu, wl, wd, s, fr, qs, qi, qg, qs_cs, qi_cs, qg_cs, qxs_ini_nzone, qxi_ini_nzone, qxg_ini_nzone = read_restart(restart_file)
            w = wu + wl + wd
            qs_lag = reshape(tile(zeros(NZONE), NSTEP + p_lag[0]), (NSTEP + p_lag[0], NZONE))
            qi_lag = reshape(tile(zeros(NZONE), NSTEP + p_lag[0]), (NSTEP + p_lag[0], NZONE))
            qg_lag = reshape(tile(zeros(NZONE), NSTEP + p_lag[0]), (NSTEP + p_lag[0], NZONE))
            qs_lag = array(qs_lag).T
            qi_lag = array(qi_lag).T
            qg_lag = array(qg_lag).T
        W_nzone = zeros(NZONE)
        WU_nzone = zeros(NZONE)
        WL_nzone = zeros(NZONE)
        WD_nzone = zeros(NZONE)
        FR_nzone = zeros(NZONE)
        S_nzone = zeros(NZONE)
        QS_nzone = zeros(NZONE)
        QI_nzone = zeros(NZONE)
        QG_nzone = zeros(NZONE)
        QS_CS_nzone = zeros(NZONE)
        QI_CS_nzone = zeros(NZONE)
        QG_CS_nzone = zeros(NZONE)
        qxs_nzone = reshape(tile(zeros(NZONE), maxdp), (maxdp, NZONE)).T
        qxi_nzone = reshape(tile(zeros(NZONE), maxdp), (maxdp, NZONE)).T
        qxg_nzone = reshape(tile(zeros(NZONE), maxdp), (maxdp, NZONE)).T
        all_et_output = []
        all_runoff_yield = []
        all_wu_output = []
        all_wl_output = []
        all_wd_output = []
        print('======================================Start Simulation=====================================')
        for itime in range(NSTEP):
            et_output = []
            runoff_yield_output = []
            print('模拟时间:', date_justification[itime])
            rs = zeros(NZONE)
            rg = zeros(NZONE)
            ri = zeros(NZONE)
            qsz = 0
            qiz = 0
            qgz = 0
            for izone in range(NZONE):
                prec_i = Prec[itime][izone]
                evapo_i = Evapo[itime][izone]
                KC_i = p_kc[izone]
                WLM_i = p_wlm[izone]
                C_i = p_c[izone]
                IMP_i = p_imp[izone]
                B_i = p_b[izone]
                WUM_i = p_wum[izone]
                WM_i = p_wm[izone]
                if itime < 1:
                    W_ini_i = w[izone]
                    WU_ini_i = wu[izone]
                    WL_ini_i = wl[izone]
                    WD_ini_i = wd[izone]
                else:
                    W_ini_i = W_nzone[izone]
                    WU_ini_i = WU_nzone[izone]
                    WL_ini_i = WL_nzone[izone]
                    WD_ini_i = WD_nzone[izone]
                ET_i, EU_i, EL_i, ED_i, W_i, WU_i, WL_i, WD_i, ND_i, PED_i, RD_i, R_Yield_i = Yield(KC_i, evapo_i, WU_ini_i, WL_ini_i, WD_ini_i, C_i, B_i, WUM_i, WLM_i, WM_i, W_ini_i, IMP_i, prec_i)
                et_output.append(ET_i)
                runoff_yield_output.append(R_Yield_i)
                W_nzone[izone] = W_i
                WU_nzone[izone] = WU_i
                WL_nzone[izone] = WL_i
                WD_nzone[izone] = WD_i
                SM_i = p_sm[izone]
                EX_i = p_ex[izone]
                KG_i = p_kg[izone]
                KI_i = p_ki[izone]
                if itime < 1:
                    FR_ini_i = fr[izone]
                    S_ini_i = s[izone]
                else:
                    FR_ini_i = FR_nzone[izone]
                    S_ini_i = S_nzone[izone]
                RS_i, RI_i, RG_i, FR_i, S_i = Divide3Source(IMP_i, SM_i, EX_i, KG_i, KI_i, PED_i, RD_i, ND_i, FR_ini_i, S_ini_i)
                FR_nzone[izone] = FR_i
                S_nzone[izone] = S_i
                CG_i = p_cg[izone]
                CI_i = p_ci[izone]
                CP_i = cp[izone]
                if itime < 1:
                    QI_ini_i = qi[izone]
                    QG_ini_i = qg[izone]
                else:
                    QI_ini_i = QI_nzone[izone]
                    QG_ini_i = QG_nzone[izone]
                QS_i, QI_i, QG_i = hsRouting(CG_i, CI_i, RS_i, RI_i, RG_i, CP_i, QI_ini_i, QG_ini_i)
                QS_nzone[izone] = QS_i
                QI_nzone[izone] = QI_i
                QG_nzone[izone] = QG_i
                CS_i = p_cs[izone]
                LAG_i = p_lag[izone]
                if itime < 1:
                    QS_CS_ini_i = qs_cs[izone]
                    QI_CS_ini_i = qi_cs[izone]
                    QG_CS_ini_i = qg_cs[izone]
                else:
                    QS_CS_ini_i = QS_CS_nzone[izone]
                    QI_CS_ini_i = QI_CS_nzone[izone]
                    QG_CS_ini_i = QG_CS_nzone[izone]
                QS_CS_i, QI_CS_i, QG_CS_i = Chrouting(CS_i, LAG_i, QS_i, QI_i, QG_i, QS_CS_ini_i, QI_CS_ini_i, QG_CS_ini_i)
                QS_CS_nzone[izone] = QS_CS_i
                QI_CS_nzone[izone] = QI_CS_i
                QG_CS_nzone[izone] = QG_CS_i
                lag = p_lag[izone]
                dp_i = p_dp[izone]
                KE_i = p_ke[izone]
                XE_i = p_xe[izone]
                qs_lag[izone][itime + lag] = QS_CS_i
                qi_lag[izone][itime + lag] = QI_CS_i
                qg_lag[izone][itime + lag] = QG_CS_i
                qss_ini_i = qs_lag[izone][itime]
                qii_ini_i = qi_lag[izone][itime]
                qgg_ini_i = qg_lag[izone][itime]
                qsig_i = qss_ini_i + qii_ini_i + qgg_ini_i
                if itime < 1:
                    qxs_ini = qxs_ini_nzone[izone]
                    qxi_ini = qxi_ini_nzone[izone]
                    qxg_ini = qxg_ini_nzone[izone]
                else:
                    qxs_ini = qxs_nzone[izone]
                    qxi_ini = qxi_nzone[izone]
                    qxg_ini = qxg_nzone[izone]
                qss_i, qxs = Muskingum(dp_i, KE_i, XE_i, DT, qss_ini_i, qxs_ini)
                qii_i, qxi = Muskingum(dp_i, KE_i, XE_i, DT, qii_ini_i, qxi_ini)
                qgg_i, qxg = Muskingum(dp_i, KE_i, XE_i, DT, qgg_ini_i, qxg_ini)
                qxs_nzone[izone] = qxs
                qxi_nzone[izone] = qxi
                qxg_nzone[izone] = qxg
                qsz = qsz + qss_i
                qiz = qiz + qii_i
                qgz = qgz + qgg_i
            all_et_output.append(et_output)
            all_runoff_yield.append(runoff_yield_output)
            all_wu_output.append(WU_nzone)
            all_wl_output.append(WL_nzone)
            all_wd_output.append(WD_nzone)
            sum_qg = qgz
            sum_qig = qgz + qiz
            sum_qsig = qgz + qiz + qsz
            all_qg.append(sum_qg)
            all_qig.append(sum_qig)
            all_qsig.append(sum_qsig)

            # Write restart file, need to be debugged
            if date_justification[itime].hour == 00 or date_justification[itime] == end_date:
                restart_date = date_justification[itime].strftime('%Y%m%d%H')
                print(f'{restart_date} is the restart time')
                write_restart(zones_columns, WU_nzone, WL_nzone, WD_nzone, S_nzone, FR_nzone, QS_nzone, QI_nzone, QG_nzone, QS_CS_nzone, QI_CS_nzone, QG_CS_nzone, qxs_nzone, qxi_nzone, qxg_nzone, restart_dir, restart_date)
            
        print('=============================Simulation has been completed=============================')
        print('Start to output ET results...')
        all_et_output = DataFrame(all_et_output)
        all_et_output.columns = zones_columns
        all_et_output.insert(0, 'Date', modeling_date_range)
        all_et_output.to_csv(Result_Dir + '/ET.csv', index=False)
        print('ET results output completed')
        print('Start to output runoff yield results...')
        all_runoff_yield = DataFrame(all_runoff_yield)
        all_runoff_yield.columns = zones_columns
        all_runoff_yield.insert(0, 'Date', modeling_date_range)
        all_runoff_yield.to_csv(Result_Dir + '/Runoff_Yield.csv', index=False)
        print('runoff yield results output completed')
        print('Start to output WU results...')
        all_wu_output = DataFrame(all_wu_output)
        all_wu_output.columns = zones_columns
        all_wu_output.insert(0, 'Date', modeling_date_range)
        all_wu_output.to_csv(Result_Dir + '/WU.csv', index=False)
        print('WU results output completed')
        print('Start to output WL results...')
        all_wl_output = DataFrame(all_wl_output)
        all_wl_output.columns = zones_columns
        all_wl_output.insert(0, 'Date', modeling_date_range)
        all_wl_output.to_csv(Result_Dir + '/WL.csv', index=False)
        print('WL results output completed')
        print('Start to output WD results...')
        all_wd_output = DataFrame(all_wd_output)
        all_wd_output.columns = zones_columns
        all_wd_output.insert(0, 'Date', modeling_date_range)
        all_wd_output.to_csv(Result_Dir + '/WD.csv', index=False)
        print('WD results output completed')
        print('Start to output streamflow results...')
        all_qsig = DataFrame(all_qsig)
        all_qsig.columns = ['streamflow']
        all_qsig.insert(0, 'Date', modeling_date_range)
        all_qsig.to_csv(Result_Dir + '/streamflow.csv', index=False)
        print('Streamflow results output completed')
        print('=============================Simulation results have been output=============================')
    except Exception as e:
        print(f'Error during simulation: {e}0')
    finally:
        input('Press any key to exit')