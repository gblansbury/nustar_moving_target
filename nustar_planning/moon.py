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
            astrometric = observer.at(t).observe(jupiter)
            ra, dec, distance = astrometric.radec()

        # Store output in degrees
            radeg = ra.to(u.deg)
            decdeg = dec.to(u.deg)


            if show is True:
                print(tstart.isoformat()+' {:.5f}  {:.5f}'.format(radeg.value, decdeg.value))

            if outfile is not None:
                f.write(tstart.isoformat()+' {:.5f}  {:.5f}'.format(radeg.value, decdeg.value)+'\n')
    
    if outfile is not None:
        f.close()
    
        
    return
    