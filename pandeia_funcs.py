'''functions for jwst/pandeia stuff'''
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import numpy as np
import copy

from dd.spectrum import bnu_wav_micron as bnu

def scene_star(sptype, mag, id=1, band='johnson,v',
               name='generic source'):
    '''Get a star to put in a scene.
        
    Parameters
    ----------
    sptype : str
        Spectral type of star.
    mag : float
        Magnitude of star in band.
    id : int
        ID number of star
    band : str, optional
        Band in which mag applies.
    name : str, optional
        Name of source.
    '''
    
    s = {}
    s['id'] = id
    s['position'] = {}
    s['position']['orientation'] = 0.0
    s['position']['x_offset'] = 0.0
    s['position']['y_offset'] = 0.0
    s['shape'] = {}
    s['shape']['geometry'] = 'point'
    s['spectrum'] = {}
    s['spectrum']['lines'] = []
    s['spectrum']['name'] = name
    s['spectrum']['normalization'] = {}
    s['spectrum']['normalization']['bandpass'] = band
    s['spectrum']['normalization']['norm_flux'] = mag
    s['spectrum']['normalization']['norm_fluxunit'] = 'vegamag'
    s['spectrum']['normalization']['type'] = 'photsys'
    s['spectrum']['redshift'] = 0.0
    s['spectrum']['sed'] = {}
    s['spectrum']['sed']['key'] = sptype
    s['spectrum']['sed']['sed_type'] = 'phoenix'

    return s


def get_dot(id=0, x=0.0, y=0.0, name='dot', norm_wave=0.0, norm_flux=0.0,
            temp=0.0, size=0.0):
    '''Get a blackbody dot to put in a scene.
    
    Parameters
    ----------
    id : int, optional
        ID number of dot.
    x : float, optional
        X position of dot.
    y : float, optional
        Y position of dot.
    name : str, optional
        Name of dot.
    norm_wave : float
        Wavelength at which normalisation applies, in micron.
    norm_flux : float
        Flux normalisation, in mJy.
    temp : float
        Temperature of blackbody.
    size : float
        Size of (circular) dot.
    '''

    s = {}
    s['id'] = id
    s['position'] = {}
    s['position']['orientation'] = 0.0
    s['position']['x_offset'] = x
    s['position']['y_offset'] = y
    s['shape'] = {}
    s['shape']['geometry'] = 'flat'
    s['shape']['major'] = size
    s['shape']['minor'] = size
    s['shape']['orientation'] = 0.0
    s['spectrum'] = {}
    s['spectrum']['lines'] = []
    s['spectrum']['name'] = name
    s['spectrum']['normalization'] = {}
    s['spectrum']['normalization']['type'] = 'at_lambda'
    s['spectrum']['normalization']['norm_wave'] = norm_wave
    s['spectrum']['normalization']['norm_flux'] = norm_flux
    s['spectrum']['normalization']['norm_waveunit'] = 'micron'
    s['spectrum']['normalization']['norm_fluxunit'] = 'mjy'
    s['spectrum']['redshift'] = 0.0
    s['spectrum']['sed'] = {}
    s['spectrum']['sed']['sed_type'] = 'blackbody'
    s['spectrum']['sed']['temp'] = temp
    return s


def add_ring(scene, r, dr, flux, temp, norm_wave, inc, pa,
             npt=None, verb=True):
    '''Add a disk to a scene, as a series of flat discs.
    
    Parameters
    ----------
    scene : list
        JWST scene dict to configure pandeia observation.
    r : float
        Radius of disk in arcsec.
    dr : float
        Width of disk in arcsec.
    flux : float
        Total flux of disk in mJy.
    temp : float
        Blackbody temperature of disk.
    norm_wave : float
        Wavelength in microns at which disk normalisation applies.
    inc : float
        Inclination of disk in degrees.
    pa : float
        Position angle of disk, E of N in degrees.
    npt : int, optional
        Number of points to make up disk.
    verb : bool, optional
        Print some useful output.
    '''

    disk = copy.deepcopy(scene[0])
    i0 = len(scene)

    disk['position']['orientation'] = 0.0
    try:
        del disk['spectrum']['sed']['key'] # unnecessary
    except:
        pass
    
    # set up a chunk of disk
    disk['shape']['geometry'] = 'flat'
    disk['shape']['major'] = dr
    disk['shape']['minor'] = dr * np.cos(np.deg2rad(inc))
    disk['shape']['orientation'] = 90 + pa

    del disk['spectrum']['normalization']['bandpass']
    disk['spectrum']['redshift'] = 0.0
    disk['spectrum']['lines'] = []
    disk['spectrum']['name'] = 'disk chunk'
    disk['spectrum']['normalization']['type'] = 'at_lambda'
    disk['spectrum']['normalization']['norm_wave'] = norm_wave
    disk['spectrum']['normalization']['norm_waveunit'] = 'micron'
    disk['spectrum']['normalization']['norm_fluxunit'] = 'mjy'

    disk['spectrum']['sed']['sed_type'] = 'blackbody'
    disk['spectrum']['sed']['temp'] = temp
    disk['spectrum']['sed']['wmin'] = 0.5
    disk['spectrum']['sed']['wmax'] = 30.0
    disk['spectrum']['sed']['sampling'] = 200.0

    if npt is None:
        npt = int(2*np.pi*r / dr)

    if verb:
        print('adding {}pt disk with r:{:g}", dr:{:g}", tot flux:{:g}mJy @ {}um, temp:{:g}K, incl:{}deg'.\
              format(npt, r, dr, flux, norm_wave, temp, inc))

    angs = np.deg2rad(np.linspace(0, 360, npt, endpoint=False))
    for i,a in enumerate(angs):
        d = copy.deepcopy(disk)
        d['id'] = i + i0+1
        d['spectrum']['normalization']['norm_flux'] = flux / float(len(angs))
        d['position']['x_offset'] = r * np.sin(a) * np.cos(np.deg2rad(inc))
        d['position']['y_offset'] = r * np.cos(a)
        scene.append(d)
        
    return scene


def add_radial_profile(scene, rad, dr, flux, temp, norm_wave, inc, pa,
                       npt=None, verb=True):
    '''Add a radial profile to a scene.
    
    Parameters
    ----------
    scene : list
        JWST scene dict to configure pandeia observation.
    rad : list
        List of annuli centers of disk in arcsec.
    dr : list
        List of widths of annuli in arcsec.
    flux : float
        List of total annulus fluxes of disk in mJy.
    temp : float
        List of blackbody temperatures of disk.
    norm_wave : float
        Wavelength in microns at which disk normalisation applies.
    inc : float
        Inclination of disk in degrees.
    pa : float
        Position angle of disk, E of N in degrees.
    npt : int, optional
        Number of points to make up disk.
    verb : bool, optional
        Print some useful output.
    '''

    for i,r in enumerate(rad):
        scene = add_ring(scene, r, dr[i], flux[i], temp[i], norm_wave,
                         inc, pa, npt=npt, verb=verb)

    return scene


def scene_spectrum(scene, first_id=0, wave=None):
    '''Get the spectrum of components in a scene.
    
    Parameters
    ----------
    scene : list
        List of scene components.
    first_id : int, optional
        ID at which to start adding up spectrum.
    wave : float or list
        Wavelengths at which to compute spectrum.
    '''

    if wave is None:
        wave = np.linspace(5,30,100)
        flux = np.zeros(wave.shape)
    else:
        try:
            flux = np.zeros(len(wave))
        except:
            flux = 0.0

    for s in scene[first_id:]:
        temp = s['spectrum']['sed']['temp']
        f = bnu(wave, temp)
        f *= s['spectrum']['normalization']['norm_flux'] / \
                bnu(s['spectrum']['normalization']['norm_wave'], temp)
        flux += f

    return wave, flux


def normalise_scene(scene, norm_flux, norm_wave=None, first_id=0):
    '''Normalise the flux of a scene, by scaling everything.
    
    Parameters
    ----------
    scene : list
        List of scene components.
    norm_flux : float
        Flux to normalise to in mJy.
    norm_wave : float, optional
        Wavelength at which normalisation applies, in micron.
    first_id : int, optional
        ID at which to start adding up spectrum.
    '''

    s = copy.deepcopy(scene)

    if norm_wave is None:
        norm_wave = s[first_id]['spectrum']['normalization']['norm_wave']

    _, tot = scene_spectrum(s, first_id=first_id, wave=norm_wave)
    norm = norm_flux / tot

    for i in range(first_id, len(scene[first_id:])+1):
        s[i]['spectrum']['normalization']['norm_flux'] *= norm

    return s


def plot_disk_scene(targ, file=None):
    '''Plot a scene with a star and a disk.
        
    Parameters
    ----------
    targ : list
        List of scene components.
    '''
    
    col = []
    stari = []
    for i,s in enumerate(targ):
        if s['shape']['geometry'] == 'flat':
            col.append( s['spectrum']['normalization']['norm_flux'] /
                        (s['shape']['major'] * s['shape']['minor']) )
        else:
            stari.append(i)

    fig,ax = plt.subplots(figsize=(5,5))
    ax.axis('equal')

    r = []
    coli = 0
    for i,s in enumerate(targ):
        if i in stari:
            ax.plot(targ[0]['position']['x_offset'], targ[0]['position']['y_offset'], '.')
        else:
            x, y = s['position']['x_offset'], s['position']['y_offset']
            r.append( np.sqrt(x*x + y*y) )
            e = Ellipse((x, y), s['shape']['major'],
                        s['shape']['minor'], angle=s['shape']['orientation'])
            ax.add_artist(e)
            e.set_facecolor( (1-np.repeat(col[coli]/np.max(col),3))/1.2 )
            coli += 1

    ax.set_ylim(-np.max(r),np.max(r))
    ax.set_xlim(-np.max(r),np.max(r))

    ax.set_xlabel('x offset / arcsec')
    ax.set_ylabel('y offset / arcsec')

    if file is not None:
        fig.savefig(file)


def show_images(ims, log=False, sub=None, title=None):
    '''Show lists of images.
    
    Parameters
    ----------
    ims : list of ndarray
        List of images (2d arrays)
    log : bool, optional
        Plot log of image.
    sub : list of ndarray, optional
         List of images to subtract from ims
    title : str, optional
        List of titles.
    '''
    plt.figure(figsize=(17,3))
    for i,im in enumerate(ims):
        plt.subplot(100 + 10*len(ims) +i+1)
        if sub is not None:
            if log:
                plt.imshow(np.log10(im - sub[i]))
            else:
                plt.imshow(im - sub[i])
        else:
            if log:
                plt.imshow(np.log10(im))
            else:                
                plt.imshow(im)
        if title is not None:
            plt.title(title[i])
        plt.colorbar()
