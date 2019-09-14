from taskinit import tb, qa, ms
import numpy as np


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
