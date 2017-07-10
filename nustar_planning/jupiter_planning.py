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

    planets = load('jup310.bsp')
    jupiter, earth = planets['jupiter'], planets['earth']
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
    
    
    return observer, jupiter, ts





def position(orbits, outfile=None,load_path=None, show=False,
    parallax_correction=False):
    '''Get the position of Jupiter at the mid-point of each orbit.
    
    Takes output of parse_occ as input.
    
    Initializes the ephemeris and then loops over each orbit to give you the
    pointing time.
    
    Currently assumed to be the midpoint of the orbit.

    Optional:

        load_path (where the bsp and SkyField data files are found.

        parllax_correction (apply the parallax correction from NuSTAR's orbit).
            Downloads the latest TLE archive from the NuSTAR SOC.

        outfile: A text file where you can store the output.
            If outfile=None then the output is written to stdout.
        
        show: Force output to stdout even if you write an output file.
        
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


    observer, jupiter, ts = init_ephem(orbits,
        load_path=load_path, show=show,
        parallax_correction=parallax_correction)

    if show is True:
        print('Aim Time            RA         Dec')

    # Loop over every orbit:
    for ind in range(len(orbits)):
        tstart = orbits.loc[ind, 'visible']
        tend = orbits.loc[ind, 'occulted']
        on_time = (tend - tstart).total_seconds()
    

        # Point at the halfway point in the orbit
        point_time = tstart + 0.5*(tend - tstart)
    
        astro_time = Time(point_time)    
        t = ts.from_astropy(astro_time)
        
        # Get the coordinates.
        astrometric = observer.at(t).observe(jupiter)
        ra, dec, distance = astrometric.radec()

        # Store output in degrees
        radeg = ra.to(u.deg)
        decdeg = dec.to(u.deg)

 #        
#         if show is True and parallax_correction is True:
#             from astropy.coordinates import SkyCoord
#             
#             radeg = ra.to(u.deg)
#             decdeg = dec.to(u.deg)
#             skyfield_ephem = SkyCoord(radeg, decdeg)
#             geocentric = earth.at(t).observe(jupiter)
#             skyfield_ephem = SkyCoord(radeg, decdeg)
#             ra2, dec2, distance2 = geocentric.radec()
#             ra2deg = ra2.to(u.deg)
#             dec2deg = dec2.to(u.deg)
# 
#             geo_ephem = SkyCoord(ra2deg, dec2deg)
#             print("Parallax corection (arcsec) {}".format(
#                 skyfield_ephem.separation(geo_ephem).arcsec))


        # Figure out how much on-target time you have:
        dt += on_time

        if show is True:
            print(tstart.isoformat()+' {:.5f}  {:.5f}'.format(radeg.value, decdeg.value))

        if outfile is not None:
            f.write(tstart.isoformat()+' {:.5f}  {:.5f}'.format(radeg.value, decdeg.value)+'\n')
    
    if outfile is not None:
        f.close()
    
    if show is True:
        print('Total accumualted time {}'.format(dt))
        
    return
    