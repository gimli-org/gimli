"""Import routines several ERT file formats."""
import re
import numpy as np
import pygimli as pg


def load(fileName, verbose=False, **kwargs):
    """Shortcut to load ERT data.

    Import Data and try to assume the file format.
    Additionally to unified data format we support the wide-spread res2dinv
    format as well as ASCII column files generated by the processing software
    of various instruments (ABEM LS, Syscal Pro, Resecs, ?)

    If this fails, install pybert and use its auto importer pybert.importData.

    Parameters
    ----------
    fileName: str

    Returns
    -------
    data: pg.DataContainer

    """
    data = pg.load(fileName)
    if isinstance(data, pg.DataContainerERT):
        return data

    try:
        pg.info("could not read unified data format for ERT ... try res2dinv")
        data = importRes2dInv(fileName)
        return data
    except Exception:
        pg.info("could not read res2dinv ... try Ascii columns")

    try:
        data = importAsciiColumns(fileName)
        return data
    except Exception as e:
        pg.info("Failed importing Ascii column file. Consider using pybert.")
        pg.info(e)

    if verbose:
        pg.info("Try to import using pybert .. if available")

    pb = pg.optImport('pybert')
    data = pb.importData(fileName)

    if isinstance(data, pg.DataContainerERT):
        return data

    pg.critical("Can't import ERT data file.", fileName)


def importRes2dInv(filename, verbose=False, return_header=False):
    """Read res2dinv format file.

    Parameters
    ----------
    filename : str
    verbose : bool [False]
    return_header : bool [False]

    Returns
    -------
    pg.DataContainerERT and (in case of return_header=True)
    header dictionary

    Format
    ------
        str - title
        float - unit spacing [m]
        int - Array Number (1-Wenner, 3-Dipole-dipole atm only)
        int - Number of Datapoints
        float - x-location given in terms of first electrode
                use 1 if mid-point location is given
        int - 0 for no IP, use 1 if IP present
        str - Phase Angle  if IP present
        str - mrad if IP present
        0,90.0 - if IP present
        dataBody
    """

    def getNonEmptyRow(i, comment='#'):
        s = next(i)
        while s[0] is comment:
            s = next(i)
        return s.split('\r\n')[0]
    # def getNonEmptyRow(...)

    with open(filename, 'r') as fi:
        content = fi.readlines()

    it = iter(content)
    header = {}
    header['name'] = getNonEmptyRow(it, comment=';')
    header['spacing'] = float(getNonEmptyRow(it, comment=';'))
    typrow = getNonEmptyRow(it, comment=';')
    typ = int(typrow.rstrip('\n').rstrip('R').rstrip('L'))

    if typ == 11:
        # independent electrode positions
        header['subtype'] = int(getNonEmptyRow(it, comment=';'))
        header['dummy'] = getNonEmptyRow(it, comment=';')
        isR = int(getNonEmptyRow(it, comment=';'))

    nData = int(getNonEmptyRow(it, comment=';'))
    xLoc = float(getNonEmptyRow(it, comment=';'))
    hasIP = int(getNonEmptyRow(it, comment=';'))

    if hasIP:
        header['ipQuantity'] = getNonEmptyRow(it, comment=';')
        header['ipUnit'] = getNonEmptyRow(it, comment=';')
        header['ipData'] = getNonEmptyRow(it, comment=';')
        ipline = header['ipData'].rstrip('\n').rstrip('\r').split(' ')
        if len(ipline) > 2:  # obviously spectral data?
            header['ipNumGates'] = int(ipline[0])
            header['ipDelay'] = float(ipline[1])
            header['onTime'] = float(ipline[-2])
            header['offTime'] = float(ipline[-1])
            header['ipDT'] = np.array(ipline[2:-2], dtype=float)
            header['ipGateT'] = np.cumsum(np.hstack((header['ipDelay'],
                                                     header['ipDT'])))

    data = pg.DataContainerERT()
    data.resize(nData)

    if typ == 9 or typ == 10:
        raise Exception("Don't know how to read:" + str(typ))

    if typ in [11, 12, 13]:  # mixed array
        res = pg.Vector(nData, 0.0)
        ip = pg.Vector(nData, 0.0)
        specIP = []

        for i in range(nData):
            vals = getNonEmptyRow(it, comment=';').replace(',', ' ').split()

            # row starts with 4
            if int(vals[0]) == 4:
                eaID = data.createSensor(pg.Pos(float(vals[1]),
                                                float(vals[2])))
                ebID = data.createSensor(pg.Pos(float(vals[3]),
                                                float(vals[4])))
                emID = data.createSensor(pg.Pos(float(vals[5]),
                                                float(vals[6])))
                enID = data.createSensor(pg.Pos(float(vals[7]),
                                                float(vals[8])))
            elif int(vals[0]) == 3:
                eaID = data.createSensor(pg.Pos(float(vals[1]),
                                                float(vals[2])))
                ebID = -1
                emID = data.createSensor(pg.Pos(float(vals[3]),
                                                float(vals[4])))
                enID = data.createSensor(pg.Pos(float(vals[5]),
                                                float(vals[6])))
            elif int(vals[0]) == 2:
                eaID = data.createSensor(pg.Pos(float(vals[1]),
                                                float(vals[2])))
                ebID = -1
                emID = data.createSensor(pg.Pos(float(vals[3]),
                                                float(vals[4])))
                enID = -1
            else:
                raise Exception('dont know how to handle row', vals[0])
            res[i] = float(vals[int(vals[0])*2+1])
            if hasIP:
                # ip[i] = float(vals[int(vals[0])*2+2])
                ipCol = int(vals[0])*2+2
                ip[i] = float(vals[ipCol])
                if 'ipNumGates' in header:
                    specIP.append(vals[ipCol:])

            data.createFourPointData(i, eaID, ebID, emID, enID)

        if isR:
            data.set('r', res)
        else:
            data.set('rhoa', res)

        if hasIP:
            data.set('ip', ip)
            if 'ipNumGates' in header:
                A = np.array(specIP, dtype=float)
                A[A > 1000] = -999
                A[A < -1000] = -999
                for i in range(header['ipNumGates']):
                    data.set('ip'+str(i+1), A[:, i])
    else:  # not type 11-13
        # amount of values per column per typ
        nntyp = [0, 3, 3, 4, 3, 3, 4, 4, 3, 0, 0, 8, 10]
        nn = nntyp[typ] + hasIP  # number of columns

        dataBody = np.zeros((nn, nData))
        for i in range(nData):
            vals = getNonEmptyRow(it, comment=';').replace(',', ' ').split()
            dataBody[:, i] = np.array(vals, dtype=float)

        XX = dataBody[0]
        EL = dataBody[1]
        SP = pg.Vector(nData, 1.0)

        if nn - hasIP == 4:
            SP = dataBody[2]

        AA = None
        BB = None
        NN = None
        MM = None

        if typ == 1:  # Wenner
            AA = XX - xLoc * EL * 1.5
            MM = AA + EL
            NN = MM + EL
            BB = NN + EL
        elif typ == 2:  # Pole-Pole
            AA = XX - xLoc * EL * 0.5
            MM = AA + EL
        elif typ == 3:  # Dipole-Dipole
            AA = XX - xLoc * EL * (SP / 2. + 1.)
            BB = AA + EL
            MM = BB + SP * EL
            NN = MM + EL
        elif typ == 3:  # Dipole-Dipole
            AA = XX - xLoc * EL * (SP / 2. + 1.)
            BB = AA + EL
            MM = BB + SP * EL
            NN = MM + EL
        elif typ == 4:  # WENNER-BETA
            AA = XX - xLoc * EL * 1.5
            BB = AA + EL
            MM = BB + EL
            NN = MM + EL
        elif typ == 5:  # WENNER-GAMMA
            AA = XX - xLoc * EL * 1.5
            MM = AA + EL
            BB = MM + EL
            NN = BB + EL
        elif typ == 6:  # POLE-DIPOLE
            AA = XX - xLoc * SP * EL - (SP - 1.) * (SP < 0.) * EL
            MM = AA + SP * EL
            NN = MM + np.sign(SP) * EL
        elif typ == 7:  # SCHLUMBERGER
            AA = XX - xLoc * EL * (SP + 0.5)
            MM = AA + SP * EL
            NN = MM + EL
            BB = NN + SP * EL
        else:
            raise Exception('Datatype ' + str(typ) + ' not yet suppoted')

        for i in range(len(AA)):
            if AA is not None:
                eaID = data.createSensor(pg.Pos(AA[i], 0.0))
            else:
                eaID = -1

            if BB is not None:
                ebID = data.createSensor(pg.Pos(BB[i], 0.0))
            else:
                ebID = -1

            if MM is not None:
                emID = data.createSensor(pg.Pos(MM[i], 0.0))
            else:
                emID = -1

            if NN is not None:
                enID = data.createSensor(pg.Pos(NN[i], 0.0))
            else:
                enID = -1

            data.createFourPointData(i, eaID, ebID, emID, enID)

        data.set('rhoa', dataBody[nn - hasIP - 1])
        if hasIP:
            data.set('ip', dataBody[nn - 1])

    istopo = int(getNonEmptyRow(it, comment=';'))
    if istopo:
        ntopo = int(getNonEmptyRow(it, comment=';'))
        ap = data.additionalPoints()
        for i in range(ntopo):
            strs = getNonEmptyRow(it, comment=';').replace(',', ' ').split()
            ap.push_back(pg.Pos([float(s) for s in strs]))

        data.setAdditionalPoints(ap)

    data.sortSensorsX()
    data.sortSensorsIndex()
    if return_header:
        return data, header
    else:
        return data
# def importRes2dInv(...)


def importAsciiColumns(filename, verbose=False, return_header=False):
    """Import any ERT data file organized in columns with column header.

    Input can be:
        * Terrameter LS or SAS Ascii Export format, e.g.
    Time MeasID DPID Channel A(x) A(y) A(z) B(x) B(y) B(z) M(x) M(y) M(z)
    N(x) N(y) N(z) F(x) F(y) F(z) Note I(mA) Uout(V) U(V)        SP(V) R(O)
    Var(%)         Rhoa Cycles Pint Pext(V) T(°C) Lat Long
    2016-09-14 07:01:56 73 7 1 8 1 1 20 1 1 12 1 1
    16 1 1 14 1 2.076  99.8757 107.892 0.0920761 0 0.921907
    0.196302 23.17 1 12.1679 12.425 42.1962 0 0
        * Resecs Output format

    """
    data = pg.DataContainerERT()
    header = {}
    with open(filename, 'r', encoding='iso-8859-15') as fi:
        content = fi.readlines()
        if content[0].startswith('Injection'):  # Resecs lead-in
            for n in range(20):
                if len(content[n]) < 2:
                    break

            content = content[n+1:]

        if content[0].startswith('Filename'):  # ABEM lead-in
            for n in range(1000):
                if content[n].find("MeasID") > 0:
                    break

            for i in range(n):
                sp = content[i].split(":")
                if len(sp) > 1:
                    tok = sp[0].lstrip("\t").lstrip("- ")
                    header[tok] = sp[1].rstrip("\n").rstrip("\r")

            for last in range(len(content)-1, -1, -1):
                if content[last].find("---") == 0:
                    print(content[last])
                    last -= 1
                    print(content[last])
                    while len(content[last]) < 3:
                        last -= 1

                    last += 1
                    break
            if last <= 1:
                last = len(content)

            content = content[n:last]

        d = readAsDictionary(content, sep='\t')
        if len(d) < 2:
            d = readAsDictionary(content)
        nData = len(next(iter(d.values())))
        data.resize(nData)
        if 'Spa.1' in d:  # Syscal Pro
            abmn = ['Spa.1', 'Spa.2', 'Spa.3', 'Spa.4']
            if verbose:
                pg.debug("detected Syscalfile format")
        elif 'A(x)' in d:  # ABEM Terrameter
            abmn = ['A', 'B', 'M', 'N']
            if verbose:
                pg.debug("detected ABEM file format")
        elif 'xA' in d:  # Workbench TX2 processed data
            abmn = ['xA', 'xB', 'xM', 'xN']
            if verbose:
                pg.debug("detected Workbench file format")
        elif 'C1(x)' in d or 'C1(xm)' in d:  # Resecs
            abmn = ['C1', 'C2', 'P1', 'P2']
            if verbose:
                pg.debug("detected RESECS file format")
        else:
            pg.debug("no electrode positions found!")
            pg.debug("Keys are:", d.keys())
            raise Exception("No electrode positions found!")
        for i in range(nData):
            if abmn[0]+'(z)' in d:
                eID = [data.createSensor([d[se+'(x)'][i], d[se+'(y)'][i],
                                          d[se+'(z)'][i]]) for se in abmn]
            elif abmn[0]+'(zm)' in d:
                eID = [data.createSensor([d[se+'(xm)'][i], d[se+'(ym)'][i],
                                          d[se+'(zm)'][i]]) for se in abmn]
            elif abmn[0]+'(y)' in d:
                eID = [data.createSensor([d[se+'(x)'][i], d[se+'(y)'][i],
                                          0.]) for se in abmn]
            elif abmn[0]+'(ym)' in d:
                eID = [data.createSensor([d[se+'(xm)'][i], d[se+'(ym)'][i],
                                          0.]) for se in abmn]
            elif abmn[0]+'(x)' in d:
                eID = [data.createSensor([d[se+'(x)'][i], 0.,
                                          0.]) for se in abmn]
            elif abmn[0]+'(xm)' in d:
                eID = [data.createSensor([d[se+'(xm)'][i], 0.,
                                          0.]) for se in abmn]
            else:
                eID = [data.createSensor([d[se][i], 0., 0.]) for se in abmn]

            data.createFourPointData(i, *eID)

        # data.save('tmp.shm', 'a b m n')
        tokenmap = {'I(mA)': 'i', 'I': 'i', 'In': 'i', 'Vp': 'u',
                    'VoltageV': 'u', 'U': 'u', 'U(V)': 'u', 'UV': 'u',
                    'R(Ohm)': 'r', 'RO': 'r', 'R(O)': 'r', 'Res': 'r',
                    'Rho': 'rhoa', 'AppROhmm': 'rhoa', 'Rho-a(Ohm-m)': 'rhoa',
                    'Rho-a(Om)': 'rhoa',
                    'Var(%)': 'err', 'D': 'err', 'Dev.': 'err', 'Dev': 'err',
                    'M': 'ma', 'P': 'ip', 'IP sum window': 'ip',
                    'Time': 't'}
        # Unit conversions (mA,mV,%), partly automatically assumed
        unitmap = {'I(mA)': 1e-3, 'Var(%)': 0.01,  # ABEM
                   'U': 1e-3, 'I': 1e-3, 'D': 0.01,  # Resecs
                   'Dev.': 0.01, 'In': 1e-3, 'Vp': 1e-3}  # Syscal
        abmn = ['a', 'b', 'm', 'n']
        if 'Cycles' in d:
            d['stacks'] = d['Cycles']
        for key in d.keys():
            vals = np.asarray(d[key])
            if key.startswith('IP sum window'):  # there is a trailing number
                key = 'IP sum window'  # apparently not working
            if np.issubdtype(vals.dtype, np.floating,  # 'float'  'int'
                             ) or np.issubdtype(vals.dtype, np.signedinteger):
                if key in tokenmap:  # use the standard (i, u, rhoa) key
                    if key not in abmn:
                        if verbose:
                            pg.debug("Setting", tokenmap[key], "from", key)
                        data.set(tokenmap[key],
                                 vals * unitmap.get(key, 1.0))
                else:  # use the original key if not XX(x) etc.
                    if not re.search('([x-z])', key) and key not in abmn:
                        data.set(key.replace(' ', '_'), d[key])

        r = data['u'] / data['i']
        if hasattr(d, 'R(0)'):
            if np.linalg.norm(r-d['R(O)']) < 1e4:  # no idea what's that for
                data.set('r', r)
            else:
                pg.debug("Warning! File inconsistent")

    data.sortSensorsX()
    if return_header:
        return data, header
    else:
        return data
# def importAsciiColumns(...)


def readAsDictionary(content, token=None, sep=None):  # obsolote due to numpy?
    """Read list of strings from a file as column separated dictionary.

        e.g.
        token1 token2 token3 token4
        va1    va2    val3   val4
        va1    va2    val3   val4
        va1    va2    val3   val4

    Parameters
    ----------
    content: [string]
        List of strings read from file:
        e.g.
        with open(filename, 'r') as fi:
            content = fi.readlines()
        fi.close()
    token: [string]
        If given the tokens will be the keys of the resulting dictionary.
        When token is None, tokens will be the first row values.
        When token is a empty list, the tokens will be autonamed to
        'col' + str(ColNumber)
    ret: dictionary
        Dictionary of all data
    """
    data = dict()

    if token is None:
        header = content[0].splitlines()[0].split(sep)
        token = []

        for i, tok in enumerate(header):
            tok = tok.lstrip()
            token.append(tok)

    for i, row in enumerate(content[1:]):
        vals = row.splitlines()[0].split(sep)
        for j, v in enumerate(vals):
            v = v.replace(',', '.')

            if len(token) < j+1:
                token.append('col' + str(j))
            if token[j] not in data:
                data[token[j]] = [None] * (len(content)-1)
            try:
                data[token[j]][i] = float(v)
            except Exception:
                if len(v) == 1 and v[0] == '-':
                    v = 0.0
                data[token[j]][i] = v

    return data
