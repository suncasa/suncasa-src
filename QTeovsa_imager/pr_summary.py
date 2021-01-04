def pr_summary(ms):
    ''' Given an ms (already open), print a concise informational block of text with key information.
    '''
    from astropy.coordinates import ICRS
    from astropy.time import Time
    from astropy import units as u
    #mssum = ms.summary(verbose=False, wantreturn=True)
    mdata = ms.metadata()
    trange = mdata.timerangeforobs(0)
    tstart = Time(trange['begin']['m0']['value'],format='mjd')  # Start time of file
    tend = Time(trange['end']['m0']['value'],format='mjd')      # End time of file
    dur = int(((tend-tstart)*86400).value+0.5)                  # Duration of entire file
    fieldnames = mdata.fieldnames()
    lines = []
    lines.append('Measurement Set Name: '+ms.name())
    lines.append(' ')
    lines.append('Timerange: '+tstart.iso+' : '+tend.iso+'   Duration: '+str(dur)+'s')
    pcenter = mdata.phasecenter()
    c = ICRS(pcenter['m0']['value']*u.rad, pcenter['m1']['value']*u.rad)
    lines.append('Source: '+fieldnames[0]+'   Coordinate RA: '+str(c.ra.to_string(u.hour))+' Dec: '+str(c.dec))
    scans = ms.getscansummary()
    i = 0
    lines.append(' ')
    lines.append('Scan    Source              Start Time                         End Time                    Int Time')
    while(1):
        try:
            scan = scans[str(i)]['0']
            lines.append(str(i)+'         '+fieldnames[scan['FieldId']]+'      '+Time(scan['BeginTime'],format='mjd').iso+' : '+Time(scan['EndTime'],format='mjd').iso+'       '+str(scan['IntegrationTime'])+'s')
            i += 1
        except:
            break
    spw = ms.getspectralwindowinfo()
    spwnames = mdata.namesforspws()
    i = 0
    lines.append(' ')
    lines.append('SpwID    Name   #subchan     Ch0(MHz)   ChanWid(MHz)    CtrFreq(MHz)')
    while(1):
        try:
            sp = spw[str(i)]
            ch0 = sp['Chan1Freq']/1e6
            nchan = sp['NumChan']
            df = sp['ChanWidth']/1e6
            ctrfreq = ch0 + (nchan-1)*df/2.
            lines.append(('{:5d}       '+spwnames[i]+'     {:7d}     {:10.4f}      {:8.4f}        {:10.4f}').format(i,nchan,ch0,df, ctrfreq))
            i += 1
        except:
            break
    return lines
