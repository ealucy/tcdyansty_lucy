import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
import sys, getopt
import metpy.calc as mpcalc
from metpy.plots import Hodograph, SkewT
from metpy.units import units

exp_parser = argparse.ArgumentParser()
exp_parser.add_argument('--input',  action='store', type=str, required=True)
exp_parser.add_argument('--output', action='store', type=str)
exp_parser.add_argument('--pressure', action='store', type=float)

args = exp_parser.parse_args()

try:
   infile = args.input
except:
   print("script usage: dry_sound_hw.py --input=name_of_file")
   sys.exit(2)

if args.output:
   outfile = args.output
else:
   outfile = 'output'

if args.pressure:
   parstart = args.pressure
else:
   parstart = 0.

print('Input File is:', infile)
print('Output Prefix is:', outfile)

col_names = ['pressure', 'height', 'temperature', 'dewpoint', 'direction', 'speed', 'theta', 'thetae']
df = pd.read_fwf(infile, skiprows=9, usecols=[0, 1, 2, 3, 6, 7, 8, 9], names=col_names)

p = df['pressure'].values * units.hPa
T = df['temperature'].values * units.degC
Td = (df['dewpoint'].values * 0.0 - 100.) * units.degC
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


prof = mpcalc.parcel_profile(p[pardex:], T[pardex], -100.0 * units.degC).to('degC')
npro = prof.size

buoy = []
for k in range(npro):
  buoy.append(9.81 * (prof[k]-T[pardex+k]) / T[pardex+k])


fig = plt.figure(figsize=(6,6))
skew = SkewT(fig, rotation=45)

skew.plot(p, T, 'r', linewidth=2)
skew.plot(p[pardex:], prof, 'k', linewidth=2)

ymax = np.amax([1000, p[0] / units.hPa])
skew.ax.set_ylim(ymax, 100.)
skew.ax.set_xlim(-50, 30)

adiab = np.arange(233.15, 503.15, 10) * units.degK

skew.plot_dry_adiabats(t0=adiab)

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

ylab = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]

plt.plot(dtdz, np.log(df['pressure'].values), 'k-', linewidth=2.0)
plt.axis([-2., 10., np.log(ymax), np.log(100)])
plt.yticks(np.log(ylab), ('100', '200', '300', '400', '500', '600', '700', '800', '900', '1000'))

plt.title('')
plt.xlabel('Lapse Rate (K/km)')
plt.ylabel('hectopascals')

plt.savefig(outfile + '_dtdz.png',format='png',dpi=150)
plt.close(fig)


n2 = mpcalc.brunt_vaisala_frequency_squared(z, th)

fig = plt.figure(figsize=(6,6))

ylab = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]

plt.plot(n2, np.log(df['pressure'].values), 'k-', linewidth=2.0)
plt.axis([0., 0.0015, np.log(ymax), np.log(100)])
plt.yticks(np.log(ylab), ('100', '200', '300', '400', '500', '600', '700', '800', '900', '1000'))

plt.title('')
plt.xlabel('Brunt-Vaisalla Frequency (1/s)')
plt.ylabel('hectopascals')

plt.savefig(outfile + '_bruntv.png',format='png',dpi=150)
plt.close(fig)


fig = plt.figure(figsize=(6,6))

ylab = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]

plt.plot(buoy, np.log(p[pardex:]/units.hPa), 'k-', linewidth=2.0)
plt.axis([-3.5, 0.2, np.log(ymax), np.log(100)])
plt.yticks(np.log(ylab), ('100', '200', '300', '400', '500', '600', '700', '800', '900', '1000'))

plt.title('')
plt.xlabel('Buoyancy Acceleration (m/s2)')
plt.ylabel('hectopascals')

plt.savefig(outfile + '_buoy.png',format='png',dpi=150)
plt.close(fig)

