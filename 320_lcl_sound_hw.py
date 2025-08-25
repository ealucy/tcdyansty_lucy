import matplotlib.pyplot as plt
#import matplotlib.ticker
from matplotlib.ticker import StrMethodFormatter
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
import pandas as pd
import argparse
import sys, getopt
import os

import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.plots import Hodograph, SkewT
from metpy.units import units


exp_parser = argparse.ArgumentParser()
exp_parser.add_argument('--input',  action='store', type=str, required=True)
exp_parser.add_argument('--output', action='store', type=str)
exp_parser.add_argument('--pressure', action='store', type=float)
exp_parser.add_argument('--temp', action='store', type=float)
exp_parser.add_argument('--dewpoint', action='store', type=float)

args = exp_parser.parse_args()

try:
   infile = args.input
except:
   print("script usage: lcl_sound_hw.py --input=name_of_file")
   sys.exit(2)

if args.output:
   outfile = args.output
else:
   outfile = 'output'

if args.pressure:
   parstart = args.pressure
else:
   parstart = 0.

if args.temp:
   tpert = args.temp * units.degK
else:
   tpert = 0.0001 * units.degK

if args.dewpoint:
   dpert = args.dewpoint * units.degK
else:   
   dpert = 0.0001 * units.degK


print('Input File is:', infile)
print('Output Prefix is:', outfile)

col_names = ['pressure', 'height', 'temperature', 'dewpoint', 'direction', 'speed', 'theta', 'thetae']

df = pd.read_fwf(infile, skiprows=9,
                 usecols=[0, 1, 2, 3, 6, 7, 8, 9], names=col_names)

p = df['pressure'].values * units.hPa
T = df['temperature'].values * units.degC
Td = df['dewpoint'].values * units.degC
th = df['theta'].values * units.degK
z  = df['height'].values * units.meters
wind_speed = df['speed'].values * units.knots
wind_dir = df['direction'].values * units.degrees
u, v = mpcalc.wind_components(wind_speed, wind_dir)

df['u_wind'], df['v_wind'] = mpcalc.wind_components(wind_speed, wind_dir)

# Drop any rows with all NaN values for T, Td, winds
df = df.dropna(subset=('temperature', 'dewpoint', 'direction', 'speed',
                       'u_wind', 'v_wind'), how='all').reset_index(drop=True)

nz = T.size

tvec = df['temperature'].values

if parstart > 0.0:
  parstart = parstart * units.hPa
  pdiff = 1000. * units.hPa
  for i in range(nz-1):
    if abs(p[i]-parstart) <= pdiff:
      pardex = i
      pdiff = abs(p[i]-parstart)
else:
  pardex = 0


dtdz = np.ones(T.shape) * units('delta_degree_Celsius / meter')

dtdz[0] = -(T[1]-T[0])/(z[1]-z[0]) * 1000.

for i in range(nz-2):
  dtdz[i+1] = -0.5*((T[i+1]-T[i])/(z[i+1]-z[i])+(T[i+2]-T[i+1])/(z[i+2]-z[i+1]))*1000.

dtdz[nz-1] = -(T[nz-1]-T[nz-2])/(z[nz-1]-z[nz-2])*1000.


cprof = mpcalc.parcel_profile(p[pardex:], T[pardex], Td[pardex]).to('degC')
npro = cprof.size

clcl_pressure, clcl_temperature = mpcalc.lcl(p[pardex], T[pardex], Td[pardex])

nzp = cprof.size
for i in range(nzp-1):
  if clcl_pressure-p[i] > 0. * units.hPa:
     cpdex = i+1
     break

Tnew = (T[pardex].to('degK') + tpert).to('degC')
Dnew = (Td[pardex].to('degK') + dpert).to('degC')

pprof = mpcalc.parcel_profile(p[pardex:], Tnew, Dnew).to('degC')
npro = pprof.size

plcl_pressure, plcl_temperature = mpcalc.lcl(p[pardex], Tnew, Dnew)

nzp = cprof.size
for i in range(nzp-1):
  if plcl_pressure-p[i] > 0. * units.hPa:
     ppdex = i+1
     break

fig = plt.figure(figsize=(6,6))
skew = SkewT(fig, rotation=15)

print('Sounding LCL is: ' + ('%5.1f' % (clcl_pressure / units.hPa)) + ' hPa')

skew.plot(p, T, 'r', linewidth=2)
skew.plot(p, Td, 'g', linewidth=2)
skew.plot(p[pardex:(pardex+cpdex)], cprof[0:cpdex], 'k-', linewidth=2)
skew.plot(clcl_pressure, clcl_temperature, 'ko', markerfacecolor='black')

if abs(plcl_pressure-clcl_pressure) > 0.01 * units.hPa:

  print('Perturbed LCL is: ' + ('%5.1f' % (plcl_pressure / units.hPa)) + ' hPa')
  skew.plot(p[pardex:(pardex+ppdex)], pprof[0:ppdex], 'm-', linewidth=2)
  skew.plot(plcl_pressure, plcl_temperature, 'mo', markerfacecolor='magenta')
  
  wpert = mpcalc.mixing_ratio_from_relative_humidity(plcl_pressure, plcl_temperature, 1.0)
  skew.plot_mixing_lines(mixing_ratio=np.array([wpert]).reshape(-1, 1), colors="magenta", linestyles='dashed')

ymax = np.amax([1000, p[0] / units.hPa])
skew.ax.set_ylim(ymax, 600.)
skew.ax.set_xlim(-30, 30)

adiab = np.arange(233.15, 503.15, 10) * units.degK

skew.plot_dry_adiabats(t0=adiab)
skew.plot_mixing_lines()

wcont = mpcalc.mixing_ratio_from_relative_humidity(clcl_pressure, clcl_temperature, 1.0)

skew.plot_mixing_lines(mixing_ratio=np.array([wcont]).reshape(-1, 1), colors="k", linestyles='dashed')

plt.savefig(outfile + '_skewt.png',format='png',dpi=150)
plt.close(fig)
