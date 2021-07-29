# CASAIRING - A task for computing/plotting radial image profiles.
# Copyright (c) Nordic ARC Node (2015).
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>,
# or write to the Free Software Foundation, Inc.,
# 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met:
#
# a. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# b. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the
#    distribution.
# c. Neither the name of the author nor the names of contributors may
#    be used to endorse or promote products derived from this software
#    without specific prior written permission.
#
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

import numpy as np
import pylab as pl
import os
import sys
from taskinit import gentools
ia = gentools(['ia'])[0]

def casairing(image='', chan0=0, nchan=-1, center=[-1, -1], rmax=-1.0, nrad=100, polchan=0,
              resultfile='', ncontour=10, errorbar=1.0, angle=[[0, 360]]):

    rad2deg = 180./np.pi
    rad2as = rad2deg*3600.

    image = str(image)

    print('\n\n   IRING FOR CASA - VERSION 1.1\n\n')
    print('Will use image %s\n' % image)

    try:
        success = ia.open(image)
        if not success:
            raise Exception('%s is not a valid image')
    except Exception:
        raise Exception('ERROR in ia tool!')

    print('Reading image data')
    cube = ia.getchunk()
    imcoords = ia.coordsys().torecord()
    RADec0 = imcoords['direction0']['crval']
    pix0RADec = imcoords['direction0']['crpix']
    deltaRADec = imcoords['direction0']['cdelt']

    unit = ia.summary()['unit']
    if 'spectral2' in imcoords.keys():
        freq0 = imcoords['spectral2']['wcs']['crval']
        pix0 = imcoords['spectral2']['wcs']['crpix']
        deltanu = imcoords['spectral2']['wcs']['cdelt']
    else:
        print('No frequency info in image coordinates.')
        freq0 = ia.summary()['refval'][-1]
        pix0 = ia.summary()['refpix'][-1]
        deltanu = 0.0

    ia.close()
    imshape = np.shape(cube)
    npix = imshape[0:2]
    if len(imshape) > 2:
        nstok = imshape[2]
    else:
        nstok = 1

    if len(imshape) > 3:
        imchan = imshape[3]
    else:
        imchan = 1
        print('This seems to be a continuum image')

    if len(imshape) == 2:
        cube = cube[:, :, np.newaxis, np.newaxis]
    elif len(imshape) == 3:
        cube = cube[:, :, :, np.newaxis]

    if nchan <= 0:
        nchan = imchan - chan0

    if chan0 >= imchan:
        raise Exception(
            'chan0 is too large. There are only %i channels and chan0 = %i' % (imchan, chan0))

    if chan0 + nchan > imchan:
        raise Exception(
            'nchan is too large. There are only %i channels and chan0 = %i' % (imchan, chan0))

    if nchan <= 0:
        nchan = imchan - chan0 + 1

    print('Selecting frequency channels from %i to %i' % (chan0, chan0+nchan-1))
    if center[0] < 0 or center[1] < 0:
        center[0] = npix[0]/2
        center[1] = npix[1]/2

    Rstop = np.min(map(abs, [npix[0] - center[0], center[0],
                             npix[1] - center[1], center[1]]))*np.abs(deltaRADec[0])*rad2as
    if rmax < 0:
        print('Setting rmax to %.2e arcsec' % Rstop)
        rmax = Rstop

    if rmax > Rstop:
        raise Exception(
            'rmax is too large! The maximum value, given your center, is %.2e arcsec' % Rstop)

    rmaxpix = int(rmax/deltaRADec[0])
    freqs = np.array([freq0 + deltanu*(p-pix0) for p in range(imchan)])
    RAs = np.array([RADec0[0] + deltaRADec[0]*(p-pix0RADec[0]) for p in range(npix[0])])
    Decs = np.array([RADec0[1] + deltaRADec[1]*(p-pix0RADec[1]) for p in range(npix[1])])
    RelRA = np.array([deltaRADec[0]*(p-pix0RADec[0]) for p in range(npix[0])])
    RelDec = np.array([deltaRADec[1]*(p-pix0RADec[1]) for p in range(npix[1])])

    ri = np.linspace(0., rmax, nrad)
    rays = [[] for r in ri]
    print('Computing matrix of distances')
    distx = np.linspace(-center[0], npix[0]-center[0],
                        npix[0])*deltaRADec[0]*rad2as
    disty = np.linspace(-center[1], npix[1]-center[1],
                        npix[1])*deltaRADec[1]*rad2as
    dx2 = distx*distx
    dy2 = disty*disty
    distmatr = np.sqrt(np.outer(np.ones(len(distx)), dy2)
                       + np.outer(dx2, np.ones(len(disty))))
    angmatr = np.arctan2(np.outer(np.ones(len(distx)), disty),
                         np.outer(distx, np.ones(len(disty))))

    print('Computing azimuthal average\n')
    allangs = []
    for angs in angle:
        if angs[0] > 180.:
            angs[0] -= 360.
        if angs[1] > 180.:
            angs[1] -= 360.
        allangs.append([a*np.pi/180. for a in angs])

    for i in range(len(ri)-1):
        pixels = np.where(np.logical_and(
            distmatr >= ri[i], distmatr < ri[i+1]))
        if len(pixels) > 0:
            if len(pixels[0]) == 1:
                rays[i] = [pixels[0], pixels[1]]
            else:
                angles = angmatr[pixels]
                angmask = np.zeros(np.shape(angles), dtype=np.bool)
                for ang in allangs:
                    if ang[0] < ang[1]:
                        newmask = np.logical_and(
                            angles >= ang[0], angles <= ang[1])
                    else:
                        newmask = np.logical_or(
                            angles >= ang[0], angles <= ang[1])
                    angmask = np.logical_or(angmask, newmask)
                rays[i] = [pixels[0][angmask], pixels[1][angmask]]

    AverData = np.zeros((nrad, nchan))
    StdData = np.zeros((nrad, nchan))
    for r, ray in enumerate(rays):
        sys.stdout.write('\r  Doing radial bin %i of %i' % (r, nrad))
        sys.stdout.flush()
        if len(ray) > 0:
            if len(ray[0]) == 1:
                print('\n\n Only one pixel for radial bin # %i' % r)
            for nu in range(nchan):
                AverData[r, nu] = np.average(cube[ray[0], ray[1], polchan, nu+chan0])
                if len(ray[0]) > 1:
                    StdData[r, nu] = np.std(
                        cube[ray[0], ray[1], polchan, nu+chan0])/np.sqrt(float(len(ray[0])))
        else:
            print('\n Radial bin #%i does not have any pixel. Try to increase the angle range(s)' % r)

    StdData[np.isnan(AverData)] = 0.0
    AverData[np.isnan(AverData)] = 0.0
    if len(resultfile) > 0:
        print('\nWriting results to file %s' % resultfile)
        with open(resultfile, 'w') as ff:
            print("#    Distance (as),  Frequency (GHz),  Average (%s),  Avg. Error (%s)" % (
                unit, unit), fil=ff)
            for r in range(nrad):
                for n in range(nchan):
                    print('%.8e  %.8e  %.8e  %.8e  ' % (
                        ri[r], freqs[n+chan0]/1.e9, AverData[r, n], StdData[r, n]), file=ff)

    print('\nPlotting')
    xaxis = freqs[chan0:chan0+nchan]/(1.e+9)
    yaxis = ri
    if nchan > 1:
        contlevs = np.linspace(np.min(AverData), np.max(AverData), ncontour)
        pl.figure(figsize=(12, 6))
        pl.contourf(xaxis, yaxis, AverData, levels=contlevs)
        pl.xlabel('Frequency (GHz)')
        pl.ylabel('Distance (as)')
        bar = pl.colorbar(orientation='vertical', format='%.2e')
        bar.set_label('Azimuthal average (%s)' % unit)
        pl.axis([xaxis[0], xaxis[-1], yaxis[0], yaxis[-1]])
        pl.show()
    else:
        pl.figure(figsize=(12, 6))
        pl.plot(yaxis, AverData[:, 0], '-k')
        pl.xlabel('Distance (as)')
        pl.ylabel('Azimuthal average (%s)' % unit)
        if errorbar > 0.0:
            pl.errorbar(yaxis, AverData[:, 0], errorbar * StdData[:, 0],
                        fmt='k', linestyle='None')
        pl.show()

    return True
