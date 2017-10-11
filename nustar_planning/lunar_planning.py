def init_ephem(orbits, load_path=None, show=False,
    parallax_correction=False):
    '''Initialize Skyfield ephemeris for Jupiter BSP file
    
    Takes output of io.parse_occ as input.
    
    Requires Astropy and SkyField
    
    Optional:
    
    load_path (where the bsp and SkyField data files are found.
    
    parllax_correction (apply the parallax correction from NuSTAR's orbit).
        Downloads the latest TLE archive from the NuSTAR SOC.
    
    
    Returns:
    observer, jupiter, ts
    
    The first two are Skyfield objects. The third is the Skyfield time series
    object.
        
    '''
    from skyfield.api import Loader, EarthSatellite
    from astropy.time import Time

    if load_path is None:
        load_path = './'
        load=Loader(load_path)
    else:
        load=Loader(load_path)

    planets = load('de436.bsp')
    moon, earth = planets['moon'], planets['earth']
    ts = load.timescale()

    if parallax_correction is False:
        observer = earth
    else:
        import nustar_planning.io as io        
        start_date = orbits.loc[0, 'visible']

        utc = Time(start_date)
        tlefile = io.download_tle(outdir=load_path)
        mindt, line1, line2 = io.get_epoch_tle(utc, tlefile)
        nustar = EarthSatellite(line1, line2)
        observer = earth + nustar
    
    
    return observer, moon, ts





def position(orbits, outfile=None,load_path=None, show=False,
    parallax_correction=False, steps=5):
    '''Get the instantaious position of the Moon at a number of intervals though the
    orbit.
    
    Takes output of parse_occ as input.
    
    Initializes the ephemeris and then loops over each orbit, splits the orbit up into a
    number of intervals (default is 5) to give you the instantaneous astrometric position
    of the Moon at each time.
    
    Optional:

        load_path (where the bsp and SkyField data files are found.

        parllax_correction (apply the parallax correction from NuSTAR's orbit).
            Downloads the latest TLE archive from the NuSTAR SOC.

        outfile: A text file where you can store the output.
            If outfile=None then the output is written to stdout.
        
        show: Force output to stdout even if you write an output file.
        
        steps: Number of intervals to use (default is 5).
        
        returns 
   
    '''
    from astropy.time import Time
    from datetime import timedelta
    import astropy.units as u
    

    if outfile is None and show is False:
        show=True
    
    dt = 0.
    if outfile is not None:
        f = open(outfile, 'w')
        f.write('Aim Time            RA        Dec\n')


    observer, moon, ts = init_ephem(orbits,
        load_path=load_path, show=show,
        parallax_correction=parallax_correction)

    if show is True:
        print('Aim Time            RA         Dec')

    # Loop over every orbit:
    for ind in range(len(orbits)):
        tstart = orbits.loc[ind, 'visible']
        tend = orbits.loc[ind, 'occulted']
        on_time = (tend - tstart).total_seconds()
    
        
        dt = ( on_time ) / steps
        for i in range(steps):
            point_time = tstart + timedelta(seconds=dt * i)

        
            astro_time = Time(point_time)    
            t = ts.from_astropy(astro_time)
        
        # Get the coordinates.
            astrometric = observer.at(t).observe(moon)
            ra, dec, distance = astrometric.radec()

        # Store output in degrees
            radeg = ra.to(u.deg)
            decdeg = dec.to(u.deg)


            if show is True:
                print(point_time.isoformat()+' {:.5f}  {:.5f}'.format(radeg.value, decdeg.value))

            if outfile is not None:
                f.write(point_time.isoformat()+' {:.5f}  {:.5f}'.format(radeg.value, decdeg.value)+'\n')
    
    if outfile is not None:
        f.close()
    
        
    return
    


def position_shift(orbits, outfile=None,load_path=None, show=True,
    parallax_correction=True, **kwargs):
    """Get the pointing position for the lunar limb observations.


    Takes output of parse_occ as input. Initializes the ephemeris and then loops over
    each orbit, splitting the orbit up into a number of dwells where you keep the limb of
    the Moon relatively fixed in the FoV.
    
    The default here is to move the *center* of the FoV so that it's sitting on the
    lunar limb. We move 1 Rmoon away from the center of the Moon in such a way that
    you always keep Dets 0 and 1 occulted by the lunar limb.
    
    Right now, just writes the results out to the stdout.

    Parameters
    ----------
    load_path: str, optional
        For the Skyfield download files
    parllax_correction: bool
        Apply the NuSTAR parallax corrections (True)
    min_dwell: float with Astropy time units
        Minimum amount of time to dwell at each position. Default is 0*u.s
    min_shift: float with Astropy angle units
        Shift when the Moons move this much. Default is 360*u.arcsec
    pad_time: float with Astropy time units
        Extends the orbit from the occultation file by this amount in both directions.
        Default in 5 minutes.
    pa: float with Astropy units
        Position angle of the NuSTAR FoV. Default is 0*u.deg (co-aligned with North).
    outfile: string
        File for ascii output
    
  
    Returns
    -------

    List of pointing positions to stdout
    
    """

    from astropy.time import Time
    from datetime import timedelta
    import astropy.units as u
    import numpy as np
    from astropy.coordinates import SkyCoord
    
    
    Rmoon=kwargs.get('Rmoon', 940*u.arcsec)
    min_dwell = kwargs.get('min_dwell', 0*u.s)
    min_shift = kwargs.get('min_shift', 360*u.arcsec)
    pa = kwargs.get('pa', 0*u.deg)
    pa = pa + np.pi * u.rad
    dt = kwargs.get('dt', 5.0*u.s)
    diag = kwargs.get('diag', False)
    pad_time = kwargs.get('pad_time', 5*u.min)
#     if outfile is None and show is False:
#         show=True
#     
    if outfile is not None:
        f = open(outfile, 'w')
        f.write('Arrive By            RA        Dec\n')


    observer, moon, ts = init_ephem(orbits,
        load_path=load_path, show=show,
        parallax_correction=parallax_correction)

#    if show is True:
#        print('Aim Time            RA         Dec')




    # Loop over every orbit:
    for ind in range(len(orbits)):

        tstart = orbits.loc[ind, 'visible'] - timedelta(minutes=pad_time.to(u.min).value)
        tend = orbits.loc[ind, 'occulted'] + timedelta(minutes=pad_time.to(u.min).value)
                
        on_time = (tend - tstart).total_seconds()
    

        steps = int(on_time / dt.to(u.s).value)
        last_point = None        
        for i in range(steps):
            point_time = tstart + timedelta(seconds=dt.to(u.s).value * i)

            astro_time = Time(point_time)    
            t = ts.from_astropy(astro_time)
        
        # Get the coordinates.
            astrometric = observer.at(t).observe(moon)
            ra, dec, distance = astrometric.radec()





        # Store output in degrees
            radeg = ra.to(u.deg)
            decdeg = dec.to(u.deg)
            this_point = SkyCoord(radeg, decdeg, unit="deg")

            if last_point is not None:
                
            
            
                dshift = this_point.separation(last_point)
                dwell = point_time - last_time
                #                print(dshift.arcsec)
                if (dshift.arcsec > min_shift.to(u.arcsec).value) & (dwell.seconds > min_dwell.to(u.s).value):
                    
                    # Aim halfway between the two positions
                    aim_time = 0.5*(point_time - last_time) + last_time
                    if diag is True:
                        print('Start of dwell: '+last_time.isoformat())
                        print('End of dwell: '+point_time.isoformat())
                        print('')
                        print('Time used to aim: '+aim_time.isoformat())
                        print('Dwell Duration (sec): {:.2f} \n  Dwell Drift (arcmin): {:.2f}'.format(dwell.seconds, dshift.arcmin))
                        print('')
                    astrometric_aim = observer.at(t).observe(moon)
                    ra_aim, dec_aim, distance = astrometric.radec()
     
                    
                    dec_point = dec_aim.to(u.deg) + Rmoon.to(u.deg) * np.cos(pa)
                    ra_point = ra_aim.to(u.deg) + Rmoon.to(u.deg) * np.sin(pa) / np.cos(dec_aim.to(u.deg))
                    
                                       
                    if show is True:
                        print(last_time.isoformat()+' RA: {:.5f}  Dec: {:.5f}'.format(ra_point.value, dec_point.value))

                    if outfile is not None:
                        f.write(last_time.isoformat()+' {:.5f}  {:.5f}'.format(ra_point.value, dec_point.value)+'\n')

                    print('')
                    # Copy over for next dwell computation
                    last_time = point_time
                    last_point = this_point

                    
            else:
                last_point = this_point
                last_time = point_time
                
            
        print('')
        
    if outfile is not None:
        f.close()
    
        
    return
    