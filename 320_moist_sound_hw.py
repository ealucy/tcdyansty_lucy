import matplotlib.pyplot as plt
from matplotlib.ticker import StrMethodFormatter
import numpy as np
import pandas as pd
import argparse
import sys, getopt
import os

import metpy.calc as mpcalc
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

# We will pull the data out of the example dataset into individual variables and
# assign units.

p = df['pressure'].values * units.hPa
T = df['temperature'].values * units.degC
Td = df['dewpoint'].values * units.degC
z  = df['height'].values * units.m
wind_speed = df['speed'].values * units.knots
wind_dir = df['direction'].values * units.degrees
u, v = mpcalc.wind_components(wind_speed, wind_dir)

df['u_wind'], df['v_wind'] = mpcalc.wind_components(wind_speed, wind_dir)

# Drop any rows with all NaN values for T, Td, winds
df = df.dropna(subset=('temperature', 'dewpoint', 'direction', 'speed',
                       'u_wind', 'v_wind'), how='all').reset_index(drop=True)

th = df['theta'].values * units.degK
te = df['thetae'].values * units.degK

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


T[pardex] = (T[pardex].to('degK') + tpert).to('degC')
Td[pardex] = (Td[pardex].to('degK') + dpert).to('degC')

prof = mpcalc.parcel_profile(p[pardex:], T[pardex], Td[pardex]).to('degC')

lcl_pressure, lcl_temperature = mpcalc.lcl(p[pardex], T[pardex], Td[pardex])

prcape, prcin = mpcalc.cape_cin(p[pardex:], T[pardex:], Td[pardex:], prof)

npro = prof.size

buoy = []
for k in range(npro):
  buoy.append(9.81 * (prof[k]-T[pardex+k]) / T[pardex+k])

print('Sounding LCL is: ' + ('%5.1f' % (lcl_pressure / units.hPa)) + ' hPa')
print('Sounding CAPE is: ' + ('%5.0f' % (prcape / units.joule * units.kilogram)).strip() + ' J/kg')
print('Sounding CIN is: ' + ('%5.0f' % (prcin  / units.joule * units.kilogram)).strip() + ' J/kg')

fig = plt.figure(figsize=(6,6))
skew = SkewT(fig, rotation=45)

skew.plot(p, T, 'r', linewidth=2)
skew.plot(p, Td, 'g', linewidth=2)
skew.plot(p[pardex:], prof, 'k', linewidth=2)

ymax = np.amax([1000, p[0] / units.hPa])
skew.ax.set_ylim(ymax, 100.)
skew.ax.set_xlim(-30, 40)

adiab = np.arange(233.15, 503.15, 10) * units.degK

skew.plot_dry_adiabats(t0=adiab)
skew.plot_moist_adiabats()
skew.plot_mixing_lines(pressure=np.linspace(300, 1000) * units.mbar)

plt.savefig(outfile + '_skewt.png',format='png',dpi=150)
plt.close(fig)


fig = plt.figure(figsize=(6,6))

ylab = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]

plt.plot(th, np.log(df['pressure'].values), 'k-', linewidth=2.0)
plt.axis([275., 420., np.log(ymax), np.log(100)])
plt.yticks(np.log(ylab), ('100', '200', '300', '400', '500', '600', '700', '800', '900', '1000'))

plt.title('Potential Temperature')
plt.xlabel('theta')
plt.ylabel('hectopascals')

plt.savefig(outfile + '_theta.png',format='png',dpi=150)
plt.close(fig)


fig = plt.figure(figsize=(6,6))

plt.plot(te, np.log(df['pressure'].values), 'k-', linewidth=2.0)
plt.axis([275., 420., np.log(ymax), np.log(100)])
plt.yticks(np.log(ylab), ('100', '200', '300', '400', '500', '600', '700', '800', '900', '1000'))

plt.title('Equivalent Potential Temperature')
plt.xlabel('theta-e')
plt.ylabel('hectopascals')

plt.savefig(outfile + '_theta-e.png',format='png',dpi=150)
plt.close(fig)


fig = plt.figure(figsize=(6,6))

ylab = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]

plt.plot(dtdz, np.log(df['pressure'].values), 'k-', linewidth=2.0)
plt.axis([-2., 10., np.log(ymax), np.log(100)])
plt.yticks(np.log(ylab), ('100', '200', '300', '400', '500', '600', '700', '800', '900', '1000'))

plt.title('')
plt.xlabel('Lapse Rate (K/km)')
plt.ylabel('hectopascals')

plt.savefig(outfile + '_dtdz.png',format='png',dpi=150)
plt.close(fig)


fig = plt.figure(figsize=(6,6))

ylab = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]

plt.plot(buoy, np.log(p[pardex:]/units.hPa), 'k-', linewidth=2.0)
plt.axis([-2.5, 1.0, np.log(ymax), np.log(100)])
plt.yticks(np.log(ylab), ('100', '200', '300', '400', '500', '600', '700', '800', '900', '1000'))

plt.title('')
plt.xlabel('Buoyancy Acceleration (m/s2)')
plt.ylabel('hectopascals')

plt.savefig(outfile + '_buoy.png',format='png',dpi=150)
plt.close(fig)