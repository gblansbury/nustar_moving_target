import pandas as pd

from datetime import datetime
from astropy.time import Time


import astropy.units as u

# 
# def parse_occ(file):
#     '''Parse the occultation file that you generated usisng the orbit_model/occ script'''
#     
#     df = pd.read_csv(file, delim_whitespace=True, header=None, skiprows=6,
#                      names = ['ingress', 'ingress_ang', 'midpoint_eng', 'midpoint_ang',
#                               'egress', 'egress_ang'])
# 
#     df['visible'] = df['egress']
#     df['occulted'] = df['egress']
# 
#     for ind in range(len(df)):
#         if ind == len(df) -1:
#             break
#         df.loc[ind,('visible')] = datetime.strptime(
#                 df.loc[ind, ('egress')],
#                 '%Y:%j:%H:%M:%S')
#         df.loc[ind,('occulted')] = datetime.strptime(
#             df.loc[ind+1, ('ingress')],
#             '%Y:%j:%H:%M:%S')
#     
#     orbits = df.loc[0:len(df)-2, ('visible', 'occulted')]
#     return orbits
    

def get_jupiter_radec(orbits, outfile=None,load_path=None, show=False,
    parallax_correction=False):
    '''Get the position of Jupiter at the specified time
    
    Takes output of parse_occ as input.
    
    Requires Astropy, and SkyField
    
    Optional:
    
    load_path (where the bsp and SkyField data files are found.
    
    outfile (where you want the output to go)
    
    show (report the output). Always True if you don't specify outfile
    
    parllax_correction (apply the parallax correction from NuSTAR's orbit).
        Requries nustar_pysolar for IO routines.
        Downloads the latest TLE archive from the NuSTAR SOC.
    
    Always reports the amount of time that you've accumulated.
    
        
    '''
    
    if outfile is None and show is False:
        show=True
        
    from skyfield.api import Loader, EarthSatellite

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
        import nustar_pysolar.io as io        
        start_date = orbits.loc[0, 'visible']

        utc = Time(start_date)
        tlefile = io.download_tle(outdir=load_path)
        mindt, line1, line2 = io.get_epoch_tle(utc, tlefile)
        nustar = EarthSatellite(line1, line2)
        observer = earth + nustar


    dt = 0.
    if outfile is not None:
        f = open(outfile, 'w')
        f.write('Arrive Time          RA                  Dec\n')



    if show is True:
        print('Aim Time            RA                  Dec')

    for ind in range(len(orbits)):
        tstart = orbits.loc[ind, 'visible']
        tend = orbits.loc[ind, 'occulted']
        on_time = (tend - tstart).total_seconds()
    
        point_time = tstart + 0.5*(tend - tstart)
    
        astro_time = Time(point_time)
    
        t = ts.from_astropy(astro_time)
        
        astrometric = observer.at(t).observe(jupiter)
        ra, dec, distance = astrometric.radec()


        if show is True and parallax_correction is True:
            from astropy.coordinates import SkyCoord
            
            radeg = ra.to(u.deg)
            decdeg = dec.to(u.deg)
            skyfield_ephem = SkyCoord(radeg, decdeg)
            geocentric = earth.at(t).observe(jupiter)
            skyfield_ephem = SkyCoord(radeg, decdeg)
            ra2, dec2, distance2 = geocentric.radec()
            ra2deg = ra2.to(u.deg)
            dec2deg = dec2.to(u.deg)

            geo_ephem = SkyCoord(ra2deg, dec2deg)
            print("Parallax corection (arcsec) {}".format(
                skyfield_ephem.separation(geo_ephem).arcsec))


        radeg = ra.to(u.deg)
        decdeg = dec.to(u.deg)

        dt += on_time

        if show is True:
            print(tstart.isoformat()+' {}  {}'.format(radeg.value, decdeg.value))

        if outfile is not None:
            f.write(tstart.isoformat()+' {}  {}'.format(radeg.value, decdeg.value)+'\n')
    
    if outfile is not None:
        f.close()
    
    print('Total accumualted time {}'.format(dt))