import pandas as pd

from datetime import datetime
from astropy.time import Time





def parse_occ(file):

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
    
    df2 = df.loc[:, ('visible', 'occulted')]
    return df2

def get_jupiter_radec(orbits, outfile=None):





    if outfile is not None:
        f = open(outfile, 'w')


    from skyfield.api import Loader
    load = Loader('../data')

    ts = load.timescale()
    planets = load('jup310.bsp')
    jupiter, earth = planets['jupiter'], planets['earth']

        

    for ind in range(len(orbits)):
        tstart = orbits.loc[ind, 'visible']
        tend = orbits.loc[ind, 'occulted']
        print()
        on_time = (tend - tstart).total_seconds()
    
        point_time = tstart + 0.5*(tend - tstart)
    
        astro_time = Time(point_time)
    
        t = ts.from_astropy(astro_time)
        
        astrometric = earth.at(t).observe(jupiter)
        ra, dec, distance = astrometric.radec()
        radeg = ra.to(u.deg)
        decdeg = dec.to(u.deg)

        dt += on_time

        print(tstart.isoformat()+' {}  {}'.format(radeg.value, decdeg.value))

        if outfile is not None:
            f.write(tstart.isoformat()+' {}  {}'.format(radeg.value, decdeg.value)+'\n')
    
    if outfile is not None:
        f.close()
    
    print('Total accumualted time {}'.format(dt))