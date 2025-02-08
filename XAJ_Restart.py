# Decompiled with PyLingual (https://pylingual.io)
# Internal filename: XAJ_Restart.py
# Bytecode version: 3.10.0rc2 (3439)
# Source timestamp: 1970-01-01 00:00:00 UTC (0)

from numpy import array, zeros, tile, reshape
from pandas import read_csv, to_datetime, DataFrame

def write_restart(nzone, wu, wl, wd, s, fr, qs, qi, qg, qs_cs, qi_cs, qg_cs, qxs, qxi, qxg, restart_dir, restart_time):
    from datetime import datetime, timedelta
    restart_time = datetime.strptime(restart_time, '%Y%m%d%H')
    restart_time = restart_time + timedelta(hours=1)
    restart_time = restart_time.strftime('%Y%m%d%H')
    Muskingum = ['qxs', 'qxi', 'qxg']
    columns_id = ['wu', 'wl', 'wd', 's', 'fr', 'qs', 'qi', 'qg', 'qs_cs', 'qi_cs', 'qg_cs']
    qxs = array(qxs).T
    qxi = array(qxi).T
    qxg = array(qxg).T
    qxsig = [qxs, qxi, qxg]
    restart_list = [wu, wl, wd, s, fr, qs, qi, qg, qs_cs, qi_cs, qg_cs]
    for imkg in range(len(Muskingum)):
        for idp in range(len(qxs)):
            columns_id.append(Muskingum[imkg] + str(idp + 1).zfill(2))
            restart_list.append(qxsig[imkg][idp])
    restart_arr = array(restart_list).T
    restart_df = DataFrame(restart_arr, columns=columns_id)
    restart_df.insert(0, 'zones_id', nzone)
    restart_df.to_csv(f'{restart_dir}restart_file_{restart_time}.csv', index=False)

def read_restart(restart_file):
    restart_df = read_csv(restart_file)
    wu = array(restart_df['wu'])
    wl = array(restart_df['wl'])
    wd = array(restart_df['wd'])
    s = array(restart_df['s'])
    fr = array(restart_df['fr'])
    qs = array(restart_df['qs'])
    qi = array(restart_df['qi'])
    qg = array(restart_df['qg'])
    qs_cs = array(restart_df['qs_cs'])
    qi_cs = array(restart_df['qi_cs'])
    qg_cs = array(restart_df['qg_cs'])
    qxs = array([list(restart_df['qxs01']), list(restart_df['qxs02']), list(restart_df['qxs03']), list(restart_df['qxs04']), list(restart_df['qxs05'])])
    qxs = qxs.T
    qxi = array([list(restart_df['qxi01']), list(restart_df['qxi02']), list(restart_df['qxi03']), list(restart_df['qxi04']), list(restart_df['qxi05'])])
    qxi = qxi.T
    qxg = array([list(restart_df['qxg01']), list(restart_df['qxg02']), list(restart_df['qxg03']), list(restart_df['qxg04']), list(restart_df['qxg05'])])
    qxg = qxg.T
    return (wu, wl, wd, s, fr, qs, qi, qg, qs_cs, qi_cs, qg_cs, qxs, qxi, qxg)

def get_restart_forcing(restart_file, Prec_Dir, Evapo_Dir):
    import re
    from datetime import datetime
    Prec = read_csv(Prec_Dir)
    Evapo = read_csv(Evapo_Dir)
    Prec['Date'] = to_datetime(Prec['Date'])
    Prec.set_index('Date', inplace=True)
    zones_columns = list(Prec.columns)
    Evapo['Date'] = to_datetime(Evapo['Date'])
    Evapo.set_index('Date', inplace=True)
    match = re.search('file_(\\d+)\\.csv', restart_file)
    if match:
        number_str = match.group(1)
        date_restart = datetime.strptime(number_str, '%Y%m%d%H')
        Prec = Prec.loc[date_restart:]
        Prec = array(Prec)
        Evapo = Evapo.loc[date_restart:]
        Evapo = array(Evapo)
    else:
        date_restart = None
    return (Prec, Evapo, zones_columns, date_restart)