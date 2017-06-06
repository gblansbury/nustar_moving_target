import pandas as pd

from datetime import datetime
from astropy.time import Time


import astropy.units as u


def parse_occ(file):
    '''Parse the occultation file that you generated usisng the orbit_model/occ script'''
    
    df = pd.read_csv(file, delim_whitespace=True, header=None, skiprows=6,
                     names = ['ingress', 'ingress_ang', 'midpoint_eng', 'midpoint_ang',
                              'egress', 'egress_ang'])

    df['visible'] = df['egress']
    df['occulted'] = df['egress']

    for ind in range(len(df)):
        if ind == len(df) -1:
            break
        df.loc[ind,('visible')] = datetime.strptime(
                df.loc[ind, ('egress')],
                '%Y:%j:%H:%M:%S')
        df.loc[ind,('occulted')] = datetime.strptime(
            df.loc[ind+1, ('ingress')],
            '%Y:%j:%H:%M:%S')
    
    orbits = df.loc[0:len(df)-2, ('visible', 'occulted')]
    return orbits
    

def get_lunar_radec(tstart, tend, outfile=None,load_path=None, show=False,
    parallax_correction=False, steps=5):
    '''Get the position of Moon during an orbit.
    
    
    Inputs:
    -------
    
    tstart (datetime object): Start of the orbit
    
    tend (datetime object): End of the orbit

    steps (5.0): Number of steps to split the orbit up into. Does it evenly right now.

    parallax_correction (False): Use the correction for orbital parallax of NuSTAR.
        Requries nustar_pysolar for IO routines.
        Downloads the latest TLE archive from the NuSTAR SOC.
    
    load_path (None): Location to store the SkyField datafiles.
    
    outfile (None): Location to store the pointing positions and times.
    
    show (None): Write output. Only used if outfile is not None.

    
    Requires Astropy, and SkyField
            
    '''
    
    if outfile is None and show is False:
        show=True
    
    
    
    from skyfield.api import Loader, EarthSatellite
    from datetime import timedelta
    

    if load_path is None:
        load_path = '../data'
        load=Loader(load_path)
    else:
        load=Loader(load_path)

    planets = load('de436.bsp')
    moon, earth = planets['moon'], planets['earth']
    ts = load.timescale()

    if parallax_correction is False:
        observer = earth
    else:
        import nustar_pysolar.io as io        

        utc = Time(tstart)
        tlefile = io.download_tle(outdir=load_path)
        mindt, line1, line2 = io.get_epoch_tle(utc, tlefile)
        nustar = EarthSatellite(line1, line2)
        observer = earth + nustar


    if outfile is not None:
        f = open(outfile, 'w')
        f.write('Arrive Time          RA                  Dec\n')


    if show is True:
        print('Aim Time            RA                  Dec')

    
    dt = ( (tend - tstart).total_seconds() ) / steps
    for i in range(steps):
        point_time = tstart + timedelta(seconds=dt * i)
    
    
        astro_time = Time(point_time)
    
        t = ts.from_astropy(astro_time)
        
        astrometric = observer.at(t).observe(moon) 
        ra, dec, distance = astrometric.radec()


        if show is True and parallax_correction is True:
            from astropy.coordinates import SkyCoord
            
            radeg = ra.to(u.deg)
            decdeg = dec.to(u.deg)
            skyfield_ephem = SkyCoord(radeg, decdeg)
            geocentric = earth.at(t).observe(moon)
            skyfield_ephem = SkyCoord(radeg, decdeg)
            ra2, dec2, distance2 = geocentric.radec()
            ra2deg = ra2.to(u.deg)
            dec2deg = dec2.to(u.deg)

            geo_ephem = SkyCoord(ra2deg, dec2deg)
            print("Parallax corection (arcsec) {}".format(
                skyfield_ephem.separation(geo_ephem).arcsec))


        radeg = ra.to(u.deg)
        decdeg = dec.to(u.deg)


        if show is True:
            print(point_time.isoformat()+' {}  {}'.format(radeg.value, decdeg.value))
            print()

        if outfile is not None:
            f.write(point_time.isoformat()+' {}  {}'.format(radeg.value,decdeg.value)+'\n')

    if outfile is not None:
        f.close()
    
    return
    

# def get_lunar_steps(tstart, tend, outfile=None,load_path=None, show=False,
#     parallax_correction=False, max_move=2):
#     '''Get the slews you have to do to prevent the moon from moving too far.
#     
#     
#     Inputs:
#     -------
#     
#     tstart (datetime object): Start of the orbit
#     
#     tend (datetime object): End of the orbit
# 
#     max_move (2.0): Maximum motion (in arcminutes)
# 
#     
#     load_path (None): Location to store the SkyField datafiles.
#     
#     outfile (None): Location to store the pointing positions and times.
#     
#     show (None): Write output to a file. Only used if outfile is not None.
#     
#     
#     Notes:
#     ----
#     
#     Requires Astropy, and SkyField
#     
#     Use the correction for orbital parallax of NuSTAR.
#         Requries nustar_pysolar for IO routines.
#         Downloads the latest TLE archive from the NuSTAR SOC.
#             
#     '''
#     
#     if outfile is None and show is False:
#         show=True
#     
#     from skyfield.api import Loader, EarthSatellite
#     from datetime import timedelta
#     
# 
#     if load_path is None:
#         load_path = '../data'
#         load=Loader(load_path)
#     else:
#         load=Loader(load_path)
# 
#     planets = load('de436.bsp')
#     moon, earth = planets['moon'], planets['earth']
#     ts = load.timescale()
# 
#     if parallax_correction is False:
#         observer = earth
#     else:
#         import nustar_pysolar.io as io        
# 
#         utc = Time(tstart)
#         tlefile = io.download_tle(outdir=load_path)
#         mindt, line1, line2 = io.get_epoch_tle(utc, tlefile)
#         nustar = EarthSatellite(line1, line2)
#         observer = earth + nustar
# 
# 
#     if outfile is not None:
#         f = open(outfile, 'w')
#         f.write('Arrive Time          RA                  Dec\n')
# 
# 
#     if show is True:
#         print('Aim Time            RA                  Dec')
# 
#     
# 
#     
#     
#     dt = ( (tend - tstart).total_seconds() ) / steps
#     for i in range(steps):
#         point_time = tstart + timedelta(seconds=dt * i)
#     
#     
#         astro_time = Time(point_time)
#     
#         t = ts.from_astropy(astro_time)
#         
#         astrometric = observer.at(t).observe(moon) `
#         ra, dec, distance = astrometric.radec()
# 
# 
#         if show is True and parallax_correction is True:
#             from astropy.coordinates import SkyCoord
#             
#             radeg = ra.to(u.deg)
#             decdeg = dec.to(u.deg)
#             skyfield_ephem = SkyCoord(radeg, decdeg)
#             geocentric = earth.at(t).observe(moon)
#             skyfield_ephem = SkyCoord(radeg, decdeg)
#             ra2, dec2, distance2 = geocentric.radec()
#             ra2deg = ra2.to(u.deg)
#             dec2deg = dec2.to(u.deg)
# 
#             geo_ephem = SkyCoord(ra2deg, dec2deg)
#             print("Parallax corection (arcsec) {}".format(
#                 skyfield_ephem.separation(geo_ephem).arcsec))
# 
# 
#         radeg = ra.to(u.deg)
#         decdeg = dec.to(u.deg)
# 
# 
#         if show is True:
#             print(point_time.isoformat()+' {}  {}'.format(radeg.value, decdeg.value))
#             print()
# 
#         if outfile is not None:
#             f.write(point_time.isoformat()+' {}  {}'.format(radeg.value,decdeg.value)+'\n')
# 
#     if outfile is not None:
#         f.close()
#     
#     return