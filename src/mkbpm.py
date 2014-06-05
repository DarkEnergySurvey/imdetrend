#!/usr/bin/env python

"""
$Id$ 

Creates bad pixel mask from biascor and flatcor 

06/05/2014: Alex Drlica-Wagner <kadrlica@fnal.gov>
  Adapted J. Marriner's mkbpm.c into python.
  First commit is a functionally identical version.

@author: Alex Drlica-Wagner <kadrlica@fnal.gov>
"""

import sys
import os
import logging
import copy
import datetime
import pyfits
import numpy as np

#generic bad pixel code
BADPIX=1
#test specific bad pixel codes
BADPIX_FLAT_MIN=1
BADPIX_FLAT_MAX=2
BADPIX_FLAT_MASK=4
BADPIX_BIAS_HOT=8
BADPIX_BIAS_WARM=16
BADPIX_BIAS_MASK=32
BADPIX_BIAS_COL=64
BADPIX_EDGE=128
BADPIX_CORR=256

SHAPE=[4096,2048]

def load(filename):
    """
    Create an array of file names from single file
    or input file list.

    Parameters:
      filename : Input FITS file or list of files

    Returns:
      files    : Array of file names
    """
    try:
        # Fits file, return as list
        pyfits.open(filename)
        return np.array([filename])
    except :
        # List of files, return list
        return np.loadtxt(filename,dtype=str)
    raise IOError("Failed to load %s"%filename)

def write(bpm,filename,ccdnum):
    """
    Write the bad column mask.

    Parameters:
      bpm      : Bad pixel mask array
      filename : Output fits file
      ccdnum   : CCD number

    Returns:
    """

    hdu=pyfits.PrimaryHDU(data=copy.copy(bpm))
    hdu.scale('int16',bzero=32768,bscale=1)
    hdu.header.update('BZERO',32768,'offset data range to that of unsigned short')
    hdu.header.update('BSCALE',1,'default scaling factor')
    hdu.header.update('CCDNUM',ccdnum,'CCD Number')
    hdu.header.update('DES_EXT','MASK','Extension type')    
    time = datetime.datetime.now().ctime()
    hdu.header.update('DESMKBPM',time,'bad pixel map created')
    comment = 'DESDM: '+' '.join(sys.argv)
    hdu.header.add_comment(comment)
    hdu.writeto(filename,clobber=True)

def badcol(filename,ccdnum):
    """
    Mask the 'funky' bad columns.
    NOTE: Careful, the bad column file is indexed at one.
    
    Parameters:
    filename : Input bad column list
    ccdnum   : CCD being examined
    Returns:
    mask     : Bad column mask array
    """

    # ADW: This should be set on the command line
    MINTHRES = -5.0;

    data = np.loadtxt(filename)
    data = data[data[:,0] == ccdnum]
    dev = np.zeros(SHAPE[1]+1)
    for ccd,col,d in data:
        logging.debug("ccd=%i col=%i dev=%f"%(ccd,col,d))
        dev[col] = d

    bcol = np.zeros(SHAPE[1]+1)
    icol,ncol = 0,0

    seq = 0
    tot = 0
    # ADW: Should go arange(1,SHAPE[1]/2 + 1)
    for i in np.arange(1,SHAPE[1]/2):
        # See if part of a previous sequence
        if seq!=0 and i<=seq+2:
	  if ((dev[i]>0.0 and (dev[i-1]!=0.0 or dev[i-2]<0.0)) \
                  or (dev[i]<0.0 and dev[i-1]<0.0) \
                  or (dev[i]==0.0 and dev[i-1]<0.0)):
              seq = i
	      tot += dev[i]
              bcol[i] = True
	      continue;

        if (seq!=0 and np.fabs(tot)/size<1.0):
            logging.info("start=%i size=%f sum=%f"%(start,size,tot))

        #See if start of a new sequence
        seq = 0
        if (dev[i]>MINTHRES): continue
        seq = i
        start = seq
        size = dev[i]
        tot = size
        bcol[i] = True

    if (seq!=0 and np.fabs(tot)/size<1.0):
        logging.info("start=%i size=%f sum=%f"%(start,size,tot))

    seq = 0
    for i in np.arange(SHAPE[1]/2+1,SHAPE[1]+1)[::-1]:
        #See if part of a previous sequence
        if (seq!=0 and i>=seq-2):
            if ((dev[i]>0.0 and (dev[i+1]!=0.0 or dev[i+2]<0.0)) or \
                    (dev[i]<0.0 and dev[i+1]<0.0) or \
                    (dev[i]==0.0 and dev[i+1]<0.0)):
                seq = i
                tot += dev[i]
                bcol[i] = True 
                continue
         
        if (seq!=0 and np.fabs(tot)/size<1.0):
            logging.info("start=%i size=%f sum=%f"%(start,size,tot))

        #See if start of a new sequence
        seq = 0;
        if (dev[i]>MINTHRES): continue;
        seq = i
        start = seq
        size = dev[i]
        tot = size
        bcol[i] = True

    if (seq!=0 and np.fabs(tot)/size<1.0):
        logging.info("Bad column start=%i size=%f sum=%f"%(start,size,tot))

    naxis1 = SHAPE[0]
    naxis2 = SHAPE[1]
    logging.info("Found %i bad columns"%bcol.sum());

    for i in np.nonzero(bcol)[0]:
        logging.debug("Masking bad column=%i"%i)        
    mask = np.zeros(SHAPE,dtype=int)
    mask[:,np.nonzero(bcol)[0]] |= BADPIX_CORR

    #return bcol
    return mask;


def mkbpm(opts):
    """
    Create the bad pixel mask. Calls the bad column masking
    internally.
    
    Parameters:
    opts   : Input options (from command line)
    Returns:
    bpm    : Bad pixel mask array
    ccdnum : CCD number
    """
    # Load flats and biases
    flatnames = load(opts.flatcor)
    biasnames = load(opts.biascor)

    logging.info("Using flatcor (value>%.2f and value<%.2f) and biascor (>%.2f) to identify bad pixels",opts.flatmax,opts.flatmin,opts.biashot)

    bpm = np.zeros(SHAPE,dtype=int) 
    flatmax_bpm = np.zeros(SHAPE,dtype=int) 
    flatmin_bpm = np.zeros(SHAPE,dtype=int) 
    biashot_bpm = np.zeros(SHAPE,dtype=int) 
    biaswarm_bpm = np.zeros(SHAPE,dtype=int) 

    logging.info("Processing %i flats"%len(flatnames))
    for i,flatname in enumerate(flatnames):
        flat = pyfits.open(flatname)
        if i == 0: 
            for j,hdu in enumerate(flat):
                if hdu.header.get('DES_EXT') == 'IMAGE':
                    ext=j
                    break
            ccdnum = flat[ext].header['CCDNUM']
        else: 
            if flat[ext].header['CCDNUM'] != ccdnum:
                raise Exception("CCD number does not match")

        image = flat[ext].data
        logging.debug("Processing flat: %s"%flatname)
        flatmax_bpm += (image > opts.flatmax)
        flatmin_bpm += (image < opts.flatmin)

    logging.info("Processing %i biases"%len(biasnames))
    for i,biasname in enumerate(biasnames):
        bias = pyfits.open(biasname)
        if bias[ext].header['CCDNUM'] != ccdnum:
            raise Exception("CCD number does not match")

        image = bias[ext].data
        logging.debug("Processing bias: %s"%biasname)

        # Check that the adjacent pixels are not also hot
        nhotbias = (image > opts.biashot)
        nhotbias[1:,:]  *= (image[1:,:]>0.2*image[:-1,:])
        nhotbias[:-1,:] *= (image[:-1,:]>0.2*image[1:,:])

        biashot_bpm += nhotbias
        biaswarm_bpm += (image > opts.biaswarm)

        # Flag masks
        # if opts.flag_mask:
        #     bpm |= flat['MASK'].data

    if opts.mask_edges:
        logging.info("Masking %i edge pixels"%opts.edgesize)
        bpm |= BADPIX_EDGE
        bpm[opts.edgesize:-opts.edgesize,opts.edgesize:-opts.edgesize] &= ~BADPIX_EDGE

    #Flag bad pixels based on majority logic
    bpm |= ((flatmin_bpm >= opts.numfmin) & (bpm&BADPIX_EDGE==0))*BADPIX_FLAT_MIN
    bpm |= (flatmax_bpm >= opts.numfmax)*BADPIX_FLAT_MAX
    bpm |= (biaswarm_bpm >= opts.numbwarm)*BADPIX_BIAS_WARM
    bpm |= (biashot_bpm >= opts.numbhot)*BADPIX_BIAS_HOT

    iy,ix = np.nonzero(bpm & BADPIX_BIAS_HOT)
    # Flag the surrounding pixels as warm
    for _ix,_iy in zip(ix,iy):
        logging.info('hotpixel x=%i y=%i',_ix,_iy)
        ymin = _iy-5 if _iy-5 > 0 else 0
        ymax = _iy+6 if _iy+6 < SHAPE[0] else SHAPE[0]
        xmin = _ix-1 if _ix-1 > 0 else 0
        xmax = _ix+2 if _ix+2 < SHAPE[1] else SHAPE[1]
        bpm[ymin:ymax,xmin:xmax] |= BADPIX_BIAS_WARM

        if (ix==_ix).sum() == 1:
            logging.info("Column=%i is correctable.",_ix)
            # Flag upstream 60 pixel cross talk and downstream
            # excess (including knob?) as correctable
            if ccdnum <= 31:
                ystart = _iy+60
                if ystart < SHAPE[0]:
                    bpm[ystart::60,_ix] &= ~BADPIX_BIAS_WARM
                    bpm[ystart::60,_ix] |= BADPIX_CORR
                bpm[0:ymin,_ix] &= ~BADPIX_BIAS_WARM
                bpm[0:ymin,_ix] |= BADPIX_CORR
            else:
                ystart = _iy-60
                if ystart >= 0:
                    bpm[ystart::-60,_ix] &= ~BADPIX_BIAS_WARM
                    bpm[ystart::-60,_ix] |= BADPIX_CORR
                # ADW: Should probably be ymax
                bpm[ymin+1:,_ix] &= ~BADPIX_BIAS_WARM
                bpm[ymin+1:,_ix] |= BADPIX_CORR
        else:
            logging.info("Column=%i is NOT correctable.",_ix)
            bpm[:,_ix] |= BADPIX_BIAS_COL

    logging.info("CCDNUM=%i",ccdnum);
    colmask = badcol("funky_column.lst",ccdnum)
    bpm |= colmask

    # Clear warm bias pixels in bad columns
    bpm &= ~(BADPIX_BIAS_WARM*(bpm&BADPIX_CORR>0))
    # Also at edges
    bpm &= ~(BADPIX_BIAS_WARM*(bpm&BADPIX_EDGE>0))
    # Set generic bad pixel bit
    if (opts.doGeneric):
        bpm = (bpm > 0)*BADPIX

    return bpm,ccdnum

def diagnostics(bpm):
    """
    Count bad pixel type and print diagnostic output.
    """
    badpix = (bpm&~BADPIX_CORR>0).sum()
    corrpix=(bpm&BADPIX_CORR>0).sum()

    logging.info("Flat low=%i (%i < %g)"%((bpm&BADPIX_FLAT_MIN>0).sum(),opts.numfmin,opts.flatmin))
    logging.info("Flat high=%i (%i > %g)"%((bpm&BADPIX_FLAT_MAX>0).sum(),opts.numfmax,opts.flatmax))
    logging.info("Flat mask=%i"%(bpm&BADPIX_FLAT_MASK>0).sum())
    logging.info("Bias hot=%i (%i > %g)"%((bpm&BADPIX_BIAS_HOT>0).sum(),opts.numbhot,opts.biashot))
    logging.info("Bias warm=%i (%i > %g)"%((bpm&BADPIX_BIAS_WARM>0).sum(),opts.numbwarm,opts.biaswarm))
    logging.info("Bias mask=%i"%(bpm&BADPIX_BIAS_MASK>0).sum())
    logging.info("Bias column=%i"%(bpm&BADPIX_BIAS_COL>0).sum())
    logging.info("Edge=%i (%i)"%((bpm&BADPIX_EDGE>0).sum(),opts.edgesize))
    logging.info("Correctable=%i"%(bpm&BADPIX_CORR>0).sum())

    badfract = (100.0*badpix)/bpm.size
    corrfract = (100.0*corrpix)/bpm.size
    logging.info("Total number of uncorrectable pixels masked = %d (%.2f%%)"%(badpix,badfract))
    logging.info("Total number of correctable pixels masked = %d (%.2f%%)"%(corrpix,corrfract))

    
if __name__ == "__main__":
    import argparse
    description = "Script for generating bad pixel masks from a set of flats, biases, and science images."
    parser = argparse.ArgumentParser(description=description,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    # ADW: Uses format from old C code, should standardize with '--'
    parser.add_argument('flatcor',type=str)
    parser.add_argument('biascor',type=str)
    parser.add_argument('outfile',type=str)
    parser.add_argument('-flatmax',default=1.5,type=float,
                        help="Threshold for hot pixels in the flats.")
    parser.add_argument('-numfmax',default=1,type=int,
                        help="Number of flats where pixel must be hot.")
    parser.add_argument('-flatmin',default=0.5,type=float,
                        help="Threshold for cold pixels in the flats.")
    parser.add_argument('-numfmin',default=1,type=int,
                        help="Number of flats where pixel must be cold.")
    parser.add_argument('-biashot',default=500.,type=float,
                        help="Threshold for hot pixels in the biases.")
    parser.add_argument('-numbhot',default=1,type=int,
                        help="Number of biases where pixel must be hot.")
    parser.add_argument('-biaswarm',default=20.,type=int,
                        help="Threshold for warm pixels in the biases.")
    parser.add_argument('-numbwarm',default=1,type=int,
                        help="Number of biases where pixel must be warm.")
    parser.add_argument('-image_compare',default=None,type=str)
    parser.add_argument('-verbose',default=2,choices=range(4),type=int,
                        help='Chatter')
    parser.add_argument('-ignorebiasimage',action='store_true')
    parser.add_argument('-mask_edges',action='store_true',
                        help="Mask the edge of the image.")
    parser.add_argument('-edgesize',default=20,type=int,
                        help="Number of edge pixels to mask.")
    parser.add_argument('-respectmask',action='store_true')
    parser.add_argument('-doGeneric',action='store_true',
                        help="Set all BPM bits to BADPIX=1")
    opts = parser.parse_args()

    logging.basicConfig(level=10*(4-opts.verbose),
                        format='MKBPM:%(levelname)s: %(message)s')

    bpm,ccd = mkbpm(opts)
    diagnostics(bpm)
    write(bpm,opts.outfile,ccd)
