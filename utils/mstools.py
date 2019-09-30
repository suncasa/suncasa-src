from taskinit import tb, qa, ms
import numpy as np
from tqdm import tqdm


def get_trange(msfile):
    from astropy.time import Time
    tb.open(msfile)
    tr = np.array([tb.getcell('TIME', 0),tb.getcell('TIME', tb.nrows() - 1)]) / 24. / 3600.
    tb.close()
    return Time(tr,format='mjd')

def msclearhistory(msfile):
    '''Clears history in the a measurement sets file

    :param msfile: string
            The name of a measurement sets file

    :return:
    '''

    tb.open(msfile + '/HISTORY', nomodify=False)
    nrows = tb.nrows()
    if nrows > 0:
        tb.removerows(range(nrows))
    tb.close()


def clearflagrow(msfile, mode='clear'):
    '''

    :param msfile:
    :param mode: FLAG_ROW operation
    default: 'clear': (default) clear the FLAG_ROW
             'list': to list existing FLAG_ROW
    :return:
    '''

    if mode == 'list':
        tb.open(msfile, nomodify=True)
        a = tb.getcol('FLAG_ROW')
        nfrows = np.sum(a)
        nrows = float(len(a))
        print('{:0d} out of {:.0f} ({:.0f}%) rows are flagged in {}'.format(nfrows, nrows, nfrows / nrows * 100,
                                                                            os.path.basename(msfile)))
    elif mode == 'clear':
        tb.open(msfile, nomodify=False)
        a = tb.getcol('FLAG_ROW')
        a[:] = False
        tb.putcol('FLAG_ROW', a)
        print('reset successfully')
    tb.close()


def splitX(vis, datacolumn2='MODEL_DATA', overwrite=True, **kwargs):
    import os
    from clearcal_cli import clearcal_cli as clearcal
    from split_cli import split_cli as split
    '''

    :param msfile:
    :param datacolumn:
    :param datacolumn2:
    :return:
    '''
    kwargs2 = kwargs.copy()
    datacolumn2 = datacolumn2.upper()

    outmsfile = kwargs['outputvis']
    if os.path.exists(outmsfile):
        if overwrite:
            os.system('rm -rf {}'.format(outmsfile))
            split(vis, **kwargs)
    else:
        split(vis, **kwargs)

    for k in ['datacolumn', 'outputvis']:
        if k in kwargs2:
            kwargs2.pop(k)

    kwargs2['outputvis'] = 'tmpms.ms'
    kwargs2['datacolumn'] = datacolumn2.replace('_DATA', '')
    if os.path.exists('tmpms.ms'):
        os.system('rm -rf tmpms.ms')
    split(vis, **kwargs2)

    tb.open('tmpms.ms')
    nrows = tb.nrows()
    data = []
    for row in tqdm(range(nrows), desc='getting {} column'.format(datacolumn2), ascii=True):
        data.append(tb.getcell('DATA', row))
    tb.close()

    clearcal(outmsfile, addmodel=True)
    tb.open(outmsfile, nomodify=False)
    for row in tqdm(range(nrows), desc='writing {} column'.format(datacolumn2), ascii=True):
        tb.putcell(datacolumn2, row, data[row])
    tb.close()

    os.system('rm -rf tmpms.ms')
    return outmsfile
