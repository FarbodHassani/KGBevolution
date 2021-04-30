#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
//////////////////////////
// Utils.py
//////////////////////////
// 
// Set of functions which allows the reading and calculating
// of gevolution data within python.
//
// Requires installation of pyccl,healpy, and scipy
//
// Author: Louis Coates
//
// Last modified: January 2021
//
//////////////////////////
"""

import pyccl as ccl
import numpy as np
import healpy as hp
import struct as st
import metadata

import scipy.integrate as integrate

from multiprocessing import Pool, RawArray

import math
import time
import itertools

P_NCDM_MASS_OMEGA = 93.14         # m_ncdm/omega_ncdm [eV].
P_T_NCDM          = 0.71611       # default value for T_ncdm
HEADER_SIZE       = 264           # Size of the header in bytes.
C_PLANCK_LAW      = 4.48147e-7    # omega_g / (T_cmb [K])^4
C_BOLTZMANN_CST   = 8.61733e-5    # Boltzmann constant [eV/K]
C_SPEED_OF_LIGHT  = 2997.92458    # speed of light [100 km/s]
C_RHO_CRIT        = 2.77459457e11 # critical density [M_sun h^2 / Mpc^3]
C_FD_NORM         = 1.80308535    # Integral[q*q/(exp(q)+1), 0, infinity]


"""
//////////////////////////
// read_header
//////////////////////////
// Description:
//   Reads header of .map file to structure. (Must be next block in file)
// 
// Arguments:
//   file            File which you wish to read header from.
//   header_size     Number of bytes in the header (Should be 264).
//
// Returns:
//   header          Struct containing header information.
// 
//////////////////////////
"""

def read_header(file,header_size = metadata.HEADER_SIZE):
    s = st.Struct('I I I I ddd d d I 196c I I')
    stream=file.read(header_size)
    header = {}
    header['Nside']       = s.unpack(stream)[0]
    header['Npix']        = s.unpack(stream)[1]
    header['precision']   = s.unpack(stream)[2]
    header['Ngrid']       = s.unpack(stream)[3]
    header['direction']   = s.unpack(stream)[4:7]
    header['distance']    = s.unpack(stream)[7]
    header['boxsize']     = s.unpack(stream)[8]
    header['Nside_ring']  = s.unpack(stream)[9]
    header['header_size'] = s.unpack(stream)[206]
    header['data_size']   = s.unpack(stream)[207]
    return header

"""
//////////////////////////
// read_start_dist
//////////////////////////
// Description:
//   Reads header of .map file to structure and then returns the start distance in Mpc. (Must be next block in file)
// 
// Arguments:
//   file            File which you wish to read header from.
//   header_size     Number of bytes in the header (Should be 264).
//
// Returns:
//   distance        Distance from the observer (in Mpc) of the first HEALpix map.
// 
//////////////////////////
"""

def read_start_dist(file,header_size = metadata.HEADER_SIZE):
    
    f=open(file,'rb')
    
    f.seek(4,1)
    
    s = st.Struct('I I I I ddd d d I 196c I I')
    stream=f.read(header_size)
    header = {}
    header['Nside']       = s.unpack(stream)[0]
    header['Npix']        = s.unpack(stream)[1]
    header['precision']   = s.unpack(stream)[2]
    header['Ngrid']       = s.unpack(stream)[3]
    header['direction']   = s.unpack(stream)[4:7]
    header['distance']    = s.unpack(stream)[7]
    header['boxsize']     = s.unpack(stream)[8]
    header['Nside_ring']  = s.unpack(stream)[9]
    header['header_size'] = s.unpack(stream)[206]
    header['data_size']   = s.unpack(stream)[207]
    
    f.close()
    return header['distance'] * header['boxsize']

"""
//////////////////////////
// read_map
//////////////////////////
// Description:
//   Reads and reorders a HEALpix map from a .map file
// 
// Arguments:
//   file            File which you wish to read header from.
//   header_size     Number of bytes in the header (Should be 264).
//   skip            Number of maps to skip before starting to read.
//
// Returns:
//   pix             HEALpix map.
//   header          Struct containing header information.
// 
//////////////////////////
"""


def read_map(file,header_size = metadata.HEADER_SIZE, skip=0):
    
    f=open(file,'rb')
    
    for i in range(skip):
        f.seek(4,1)
        hdr = read_header(f,header_size)
        f.seek(hdr['data_size'],1)
        f.seek(4,1)

    f.seek(4,1)
    hdr = read_header(f,header_size)
    pixbatch_delim = np.array([0,0,0])
    pixbatch_size  = np.array([0,0,0])
    
    if (hdr['Nside_ring'] > 0 & hdr['Nside_ring'] < hdr['Nside']):
    
        #print('pixels will be reordered')
        
        if   (hdr['Npix'] <= 2 * hdr['Nside'] * (hdr['Nside'] + 1)):
            ring = math.floor((math.sqrt(2 * hdr['Npix'] + 1.01) - 1) / 2.)
        elif (hdr['Npix'] <= 2 * hdr['Nside'] * (hdr['Nside'] + 1) + 4 * (2 * hdr['Nside'] - 1) * hdr['Nside']):
            ring = ((hdr['Npix'] - 2 * hdr['Nside'] * (hdr['Nside'] + 1)) / 4 / hdr['Nside']) + hdr['Nside']
        elif (hdr['Npix'] < 12 * hdr['Nside'] * hdr['Nside']):
            ring = 12 * hdr['Nside'] * hdr['Nside'] - hdr['Npix']
            ring = math.floor((math.sqrt(2 * ring + 1.01) - 1) / 2)
            ring = 4 * hdr['Nside'] - 1 - ring
        else:
            ring = 4 * hdr['Nside'] - 1
        
        pixbatch_size[0] = (hdr['Nside'] / hdr['Nside_ring'])
        
        pixbatch_delim[1] = ring / pixbatch_size[0]
        pixbatch_delim[0] = pixbatch_delim[1]-1 if (pixbatch_delim[1] > 0) else 0
        pixbatch_delim[2] = pixbatch_delim[1]+1
        pixbatch_size[1] = (pixbatch_size[0] * (pixbatch_size[0]+1) + (2*pixbatch_size[0] - 1 - ring%pixbatch_size[0]) * (ring%pixbatch_size[0])) / 2;
        pixbatch_size[2] = ((ring%pixbatch_size[0] + 1) * (ring%pixbatch_size[0])) / 2
        pixbatch_size[0] *= pixbatch_size[0]        
        
        for p in range(3):
            if (pixbatch_delim[p] <= hdr['Nside_ring']):
                pixbatch_delim[p] = 2 * pixbatch_delim[p] * (pixbatch_delim[p]+1)
            elif (pixbatch_delim[p] <= 3 * hdr['Nside_ring']):
                pixbatch_delim[p] = 2 * hdr['Nside_ring'] * (hdr['Nside_ring']+1) + (pixbatch_delim[p]-hdr['Nside_ring']) * 4 * hdr['Nside_ring']
            elif (pixbatch_delim[p] < 4 * hdr['Nside_ring']):
                pixbatch_delim[p] = 12 * hdr['Nside_ring'] * hdr['Nside_ring'] - 2 * (4 * hdr['Nside_ring'] - 1 - pixbatch_delim[p]) * (4 * hdr['Nside_ring'] - pixbatch_delim[p])
            else:
                pixbatch_delim[p] = 12 * hdr['Nside_ring'] * hdr['Nside_ring']
    
    if   hdr['precision']==4:
        s = st.Struct(str(hdr['Npix'])+'f')
    elif hdr['precision']==8:
        s = st.Struct(str(hdr['Npix'])+'d')
        
    stream=f.read(hdr['data_size'])
    rpix = np.array(s.unpack(stream),dtype='float')
    pix=np.full(hdr['Nside']*hdr['Nside']*12,-1.6375e30,dtype='float')
    
    if (pixbatch_size[0] > 1):
        #rpix = pix
        for p in range(pixbatch_delim[0]):
            for i in range(pixbatch_size[0]):
                j = hp.ring2nest(hdr['Nside_ring'], p)
                j = j*pixbatch_size[0] + i
                j = hp.nest2ring(hdr['Nside'], j)
                pix[j] = rpix[pixbatch_size[0]*p + i]
        
        for p in range(pixbatch_delim[0],pixbatch_delim[1],1):
            q = 0
            for i in range(pixbatch_size[0]):
                j = hp.ring2nest(hdr['Nside_ring'], p)
                j = j*pixbatch_size[0] + i
                j = hp.nest2ring(hdr['Nside'], j)
                if (j < hdr['Npix']):
                    pix[j] = rpix[pixbatch_size[0]*pixbatch_delim[0] + pixbatch_size[1]*(p-pixbatch_delim[0]) + q]
                    q+=1

        for p in range(pixbatch_delim[1],pixbatch_delim[2],1):
            q = 0;
            for i in range(pixbatch_size[0]):
                j = hp.ring2nest(hdr['Nside_ring'], p)
                j = j*pixbatch_size[0] + i
                j = hp.nest2ring(hdr['Nside'], j)
                if (j < hdr['Npix']):
                    pix[j] = rpix[pixbatch_size[0]*pixbatch_delim[0] + pixbatch_size[1]*(pixbatch_delim[1]-pixbatch_delim[0]) + pixbatch_size[2]*(p-pixbatch_delim[1]) + q]
                    q+=1
    else:
        pix=rpix

    f.close()
    return pix,hdr

var_dict = {}

def init_worker(Nside,Nside_ring,pixbatch,rpix,pix):

    var_dict['Nside'] = Nside
    var_dict['Nside_ring'] = Nside_ring
    var_dict['pixbatch'] = pixbatch
    var_dict['rpix'] = rpix
    var_dict['pix'] = pix

def worker_func(params):
    p , i = params
    pix_tmp = var_dict['pix']
    rpix_tmp = var_dict['rpix']
    j = hp.ring2nest(var_dict['Nside_ring'], p)
    j = j*var_dict['pixbatch'] + i
    j = hp.nest2ring(var_dict['Nside'], j)
    pix_tmp[j] = rpix_tmp[var_dict['pixbatch']*p + i]
    
"""
//////////////////////////
// read_map_parallel
//////////////////////////
// Description:
//   Reads and reorders a HEALpix map from a .map file
// 
// Arguments:
//   file            File which you wish to read header from.
//   header_size     Number of bytes in the header (Should be 264).
//   skip            Number of maps to skip before starting to read.
//   cpus            Number of cpus available
//
// Returns:
//   pix             HEALpix map.
//   header          Struct containing header information.
// 
//////////////////////////
"""

def read_map_parallel(file,header_size= metadata.HEADER_SIZE,skip=0,cpus=4):
    
    f=open(file,'rb')
    
    for i in range(skip):
        f.seek(4,1)
        hdr = read_header(f,header_size)
        f.seek(hdr['data_size'],1)
        f.seek(4,1)

    f.seek(4,1)
    hdr = read_header(f,header_size)
    pixbatch_delim = np.array([0,0,0])
    pixbatch_size  = np.array([0,0,0])
    
    if (hdr['Nside_ring'] > 0 & hdr['Nside_ring'] < hdr['Nside']):
        
        if   (hdr['Npix'] <= 2 * hdr['Nside'] * (hdr['Nside'] + 1)):
            ring = math.floor((math.sqrt(2 * hdr['Npix'] + 1.01) - 1) / 2.)
        elif (hdr['Npix'] <= 2 * hdr['Nside'] * (hdr['Nside'] + 1) + 4 * (2 * hdr['Nside'] - 1) * hdr['Nside']):
            ring = ((hdr['Npix'] - 2 * hdr['Nside'] * (hdr['Nside'] + 1)) / 4 / hdr['Nside']) + hdr['Nside']
        elif (hdr['Npix'] < 12 * hdr['Nside'] * hdr['Nside']):
            ring = 12 * hdr['Nside'] * hdr['Nside'] - hdr['Npix']
            ring = math.floor((math.sqrt(2 * ring + 1.01) - 1) / 2)
            ring = 4 * hdr['Nside'] - 1 - ring
        else:
            ring = 4 * hdr['Nside'] - 1
        
        pixbatch_size[0] = (hdr['Nside'] / hdr['Nside_ring'])
        
        pixbatch_delim[1] = ring / pixbatch_size[0]
        pixbatch_delim[0] = pixbatch_delim[1]-1 if (pixbatch_delim[1] > 0) else 0
        pixbatch_delim[2] = pixbatch_delim[1]+1
        pixbatch_size[1] = (pixbatch_size[0] * (pixbatch_size[0]+1) + (2*pixbatch_size[0] - 1 - ring%pixbatch_size[0]) * (ring%pixbatch_size[0])) / 2;
        pixbatch_size[2] = ((ring%pixbatch_size[0] + 1) * (ring%pixbatch_size[0])) / 2
        pixbatch_size[0] *= pixbatch_size[0]        
        
        for p in range(3):
            if (pixbatch_delim[p] <= hdr['Nside_ring']):
                pixbatch_delim[p] = 2 * pixbatch_delim[p] * (pixbatch_delim[p]+1)
            elif (pixbatch_delim[p] <= 3 * hdr['Nside_ring']):
                pixbatch_delim[p] = 2 * hdr['Nside_ring'] * (hdr['Nside_ring']+1) + (pixbatch_delim[p]-hdr['Nside_ring']) * 4 * hdr['Nside_ring']
            elif (pixbatch_delim[p] < 4 * hdr['Nside_ring']):
                pixbatch_delim[p] = 12 * hdr['Nside_ring'] * hdr['Nside_ring'] - 2 * (4 * hdr['Nside_ring'] - 1 - pixbatch_delim[p]) * (4 * hdr['Nside_ring'] - pixbatch_delim[p])
            else:
                pixbatch_delim[p] = 12 * hdr['Nside_ring'] * hdr['Nside_ring']
    
    if   hdr['precision']==4:
        s = st.Struct(str(hdr['Npix'])+'f')
    elif hdr['precision']==8:
        s = st.Struct(str(hdr['Npix'])+'d')
     
    stream=f.read(hdr['data_size'])
    rpix = np.array(s.unpack(stream),dtype='float')

    if (pixbatch_size[0] > 1):

        pix = np.ctypeslib.as_ctypes(np.full(hdr['Nside']*hdr['Nside']*12,-1.6375e30,dtype='float'))
        pix_share  = RawArray('f',pix)

        param1 = range(pixbatch_delim[0])
        param2 = range(pixbatch_size[0])
        paramlist = list(itertools.product(param1,param2))

        pool=Pool(cpus,init_worker(hdr['Nside'],hdr['Nside_ring'],pixbatch_size[0],rpix,pix_share),())

        pool.map(worker_func,paramlist)
        
        pix = np.ctypeslib.as_array(pix_share)
        pool.close()

        for p in range(pixbatch_delim[0],pixbatch_delim[1],1):
            q = 0
            for i in range(pixbatch_size[0]):
                j = hp.ring2nest(hdr['Nside_ring'], p)
                j = j*pixbatch_size[0] + i
                j = hp.nest2ring(hdr['Nside'], j)
                if (j < hdr['Npix']):
                    pix[j] = rpix[pixbatch_size[0]*pixbatch_delim[0] + pixbatch_size[1]*(p-pixbatch_delim[0]) + q]
                    q+=1

        for p in range(pixbatch_delim[1],pixbatch_delim[2],1):
            q = 0;
            for i in range(pixbatch_size[0]):
                j = hp.ring2nest(hdr['Nside_ring'], p)
                j = j*pixbatch_size[0] + i
                j = hp.nest2ring(hdr['Nside'], j)
                if (j < hdr['Npix']):
                    pix[j] = rpix[pixbatch_size[0]*pixbatch_delim[0] + pixbatch_size[1]*(pixbatch_delim[1]-pixbatch_delim[0]) + pixbatch_size[2]*(p-pixbatch_delim[1]) + q]
                    q+=1
    else:
        pix=rpix

    f.close()
    return pix,hdr

"""
//////////////////////////
// get_metadata
//////////////////////////
// Description:
//   retrieves metadata from file
// 
// Arguments:
//   path          path to the file
//
// Returns: strings with path information
// 
//////////////////////////
"""


def get_metadata(path):
    lines=[]
    with open(path,'rt') as settingsfile:
        for line in settingsfile:
            lines.append(line.rstrip('\n').replace(" ","").split("=",1))
    for element in lines:
        if element[0]=="outputpath":
            output_path     = str(element[1])
        if element[0]=="genericfilebase":
            generic_file    = str(element[1])
        if element[0]=="lightconefilebase":
            lightcone_file  = str(element[1])
            
    return output_path,generic_file,lightcone_file

"""
//////////////////////////
// set_cosmology
//////////////////////////
// Description:
//   retrieves cosmology information from gevolution settings file (When there are 2 non-relativistic species)
// 
// Arguments:
//   path          path to the file
//
// Returns: Series of cosmology information within dictionaries and floats
// 
//////////////////////////
"""

def set_cosmology(path):
    lines=[]
    with open(path,'rt') as settingsfile:
        for line in settingsfile:
            lines.append(line.rstrip('\n').replace(" ","").split("=",1))

    w0_fld=-1
    Omega_fld=0

    for element in lines:
        if element[0]=="A_s":
            A_s         = float(element[1])
        if element[0]=="n_s":
            n_s         = float(element[1])
        if element[0]=="h":
            h           = float(element[1])
        if element[0]=="omega_b":
            omega_b     = float(element[1])
        if element[0]=="omega_cdm":
            omega_cdm   = float(element[1])
        if element[0]=="T_cmb":
            T_cmb       = float(element[1])
        if element[0]=="N_ur":
            N_ur        = float(element[1])
        if element[0]=="m_ncdm":
            m_ncdm      = [float(element[1].split(",",1)[0]),float(element[1].split(",",1)[1])]
        if element[0]=="w0_fld":
            w0_fld      = float(element[1])            
        if element[0]=="Omega_fld":
            Omega_fld   = float(element[1])

    Omega_g = T_cmb
    Omega_g = Omega_g * Omega_g / h
    Omega_g = Omega_g * Omega_g * C_PLANCK_LAW
    
    Omega_ur = N_ur
    Omega_ur *= (7/8) * pow(4/11, 4/3) * Omega_g
    
    Omega_rad = Omega_g + Omega_ur
    
    cosmo = ccl.Cosmology(Omega_b=(omega_b/h**2), Omega_c=(omega_cdm/h**2),h=h,A_s=A_s,n_s=n_s,T_CMB=T_cmb, Omega_g=Omega_g, w0=w0_fld)
    return omega_b,omega_cdm , h , m_ncdm , Omega_rad,Omega_fld,w0_fld, cosmo

"""
//////////////////////////
// set_cosmology_2mnu
//////////////////////////
// Description:
//   retrieves cosmology information from gevolution settings file (When there are 3 non-relativistic species)
// 
// Arguments:
//   path          path to the file
//
// Returns: Series of cosmology information within dictionaries and floats
// 
//////////////////////////
"""

def set_cosmology_3mnu(path):
    lines=[]
    with open(path,'rt') as settingsfile:
        for line in settingsfile:
            lines.append(line.rstrip('\n').replace(" ","").split("=",1))

    w0_fld=-1
    Omega_fld=0

    for element in lines:
        if element[0]=="A_s":
            A_s         = float(element[1])
        if element[0]=="n_s":
            n_s         = float(element[1])
        if element[0]=="h":
            h           = float(element[1])
        if element[0]=="omega_b":
            omega_b     = float(element[1])
        if element[0]=="omega_cdm":
            omega_cdm   = float(element[1])
        if element[0]=="T_cmb":
            T_cmb       = float(element[1])
        if element[0]=="N_ur":
            N_ur        = float(element[1])
        if element[0]=="m_ncdm":
            m_ncdm      = [float(element[1].split(",",2)[0]),float(element[1].split(",",2)[1]),float(element[1].split(",",2)[2])]


    Omega_g = T_cmb
    Omega_g = Omega_g * Omega_g / h
    Omega_g = Omega_g * Omega_g * C_PLANCK_LAW

    Omega_ur = N_ur
    Omega_ur *= (7/8) * pow(4/11, 4/3) * Omega_g

    Omega_rad = Omega_g + Omega_ur

    cosmo = ccl.Cosmology(Omega_b=(omega_b/h**2), Omega_c=(omega_cdm/h**2),h=h,A_s=A_s,n_s=n_s,T_CMB=T_cmb, Omega_g=Omega_g)
    return omega_b,omega_cdm , h , m_ncdm , Omega_rad,Omega_fld,w0_fld, cosmo

def set_lightcones(row,value_set):
    row['lightcone_num']=value_set


"""
//////////////////////////
// FermiDiracIntegral
//////////////////////////
// Description:
//   computes the integral of the relativistic Fermi-Dirac distribution
// 
// Arguments:
//   w          parameter in the F-D distribution, "(m a / kB T)^2"
//
// Returns: value for the integral
// 
//////////////////////////
"""

def FermiDiracIntegral(w):
    result = integrate.quad(lambda q: (q * q * np.sqrt(q * q + w) / (np.exp(q)+1)),0,24)
    return result[0]

"""
//////////////////////////
// bg_ncdmeach
//////////////////////////
// Description:
//   computes the background model for one ncdm species by integrating the relativistic
//   Fermi-Dirac distribution
// 
// Arguments:
//   a          scale factor at which to compute the background model
//   cosmo      structure containing the cosmological parameters
//   p          index of the ncdm species
//   m_ncdm     array of ncdm particle masses (in eV)
//
// Returns: value for the background model
// 
//////////////////////////
"""

def bg_ncdmeach(a,cosmo,p,m_ncdm):
        w = a * m_ncdm[p] / (pow(cosmo['Omega_g']* cosmo['h'] * cosmo['h'] / C_PLANCK_LAW, 0.25) * P_T_NCDM * C_BOLTZMANN_CST)
        w *= w
        return FermiDiracIntegral(w) * (m_ncdm[p] / P_NCDM_MASS_OMEGA / cosmo['h'] / cosmo['h']) * pow(cosmo['Omega_g'] * cosmo['h'] * cosmo['h'] / C_PLANCK_LAW, 0.25) * P_T_NCDM * C_BOLTZMANN_CST / m_ncdm[p] / C_FD_NORM / a

 """
//////////////////////////
// bg_ncdm
//////////////////////////
// Description:
//   computes the background model for all ncdm species by integrating the relativistic
//   Fermi-Dirac distribution
// 
// Arguments:
//   a          scale factor at which to compute the background model
//   cosmo      structure containing the cosmological parameters
//   m_ncdm     array of ncdm particle masses (in eV)
//
// Returns: value for the background model
// 
//////////////////////////
"""

def bg_ncdm(a,cosmo,m_ncdm):
    result = -1.0
    a_prev = -1.0
    if (a != a_prev):
        result = 0.0;
        a_prev = a;
        for p in range(0, len(m_ncdm)):
            result += bg_ncdmeach(a, cosmo, p, m_ncdm)
    
    return result

"""
//////////////////////////
// Hconf
//////////////////////////
// Description:
//   computes the conformal Hubble rate at given scale factor
// 
// Arguments:
//   a          scale factor
//   fourpiG    "4 pi G"
//   cosmo      structure containing the cosmological parameters
//   m_ncdm     array of ncdm particle masses (in eV)
//   Omega_rad  Dimensionful density parameter for radiation
//   Omega_fld  Dimensionful density parameter for fluid
//   w0         w0 for fluid
//
// Returns: conformal Hubble rate
// 
//////////////////////////
"""

def Hconf(a,fourpiG,cosmo,m_ncdm,Omega_rad,Omega_fld,w0_fld):
    return np.sqrt((2 * fourpiG / 3) * (((cosmo['Omega_c'] + cosmo['Omega_b'] + bg_ncdm(a, cosmo, m_ncdm)) / a) + ((1-(cosmo['Omega_c'] + cosmo['Omega_b'] + bg_ncdm(a, cosmo, m_ncdm))-Omega_rad-Omega_fld) * a * a) + (Omega_rad / a / a) + (Omega_fld/pow(a, 3 + (3 * w0_fld)))))

"""
//////////////////////////
// particleHorizon
//////////////////////////
// Description:
//   computes the particle horizon (tau) at given scale factor
// 
// Arguments:
//   a          scale factor
//   fourpiG    "4 pi G"
//   cosmo      structure containing the cosmological parameters
//   m_ncdm     array of ncdm particle masses (in eV)
//   Omega_rad  Dimensionful density parameter for radiation
//   Omega_fld  Dimensionful density parameter for fluid
//   w0         w0 for fluid
//
// Returns: particle horizon (tau)
// 
//////////////////////////
"""

def particleHorizon(a,fourpiG,cosmo,m_ncdm,Omega_rad,Omega_fld,w0_fld):
    result = integrate.quad(lambda sqrta: (2 / (sqrta * Hconf(sqrta*sqrta, 1,cosmo,m_ncdm,Omega_rad,Omega_fld,w0_fld))),np.sqrt(a) * 1.0e-7,np.sqrt(a))
    return result[0]/ np.sqrt(fourpiG)

