from taskinit import tb, qa, ms

def getAntennaPosition(vis):
    tb.open(vis + '/ANTENNA')
    position = tb.getcol('POSITION')
    diameter = tb.getcol('DISH_DIAMETER')
    antenna = tb.getcol('NAME')
    tb.close()
    return position, diameter, antenna


def getObservatoryName(ms):
    """
    Returns the observatory name in the specified ms, using the tb tool.
    -- Todd Hunter
    """
    antTable = ms + '/OBSERVATION'
    try:
        tb.open(antTable)
        myName = tb.getcell('TELESCOPE_NAME')
        tb.close()
    except:
        print("Could not open OBSERVATION table to get the telescope name: {}".format(antTable))
        myName = ''
    return (myName)


def buildConfigurationFile(vis='', cfgfile=None):
    observatory = getObservatoryName(vis)
    position, diameter, antenna = getAntennaPosition(vis)
    if cfgfile == '':
        cfgfile = vis + '.cfg'
    cfg = open(cfgfile, 'w')
    coordsys = '# coordsys=XYZ\n'
    cfg.write('# observatory={}\n'.format(observatory))
    cfg.write(coordsys)
    for s in range(len(antenna)):
        x, y, z = position[:, s]
        line = '{:.6f} {:.6f} {:.6f} {:5.1f}  {}\n'.format(x, y, z, diameter[s], antenna[s])
        cfg.write(line)
    cfg.close()
    return cfgfile
