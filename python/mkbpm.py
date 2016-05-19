#!/usr/bin/env python
"""
$Id: mkbpm.py 26607 2014-10-15 01:56:28Z kadrlica $ 

Creates bad pixel mask from biascor and flatcor 

06/05/2014: Alex Drlica-Wagner <kadrlica@fnal.gov>
  Adapted J. Marriner's mkbpm.c into python.
  First commit is a functionally identical version.
10/15/2014: Alex Drlica-Wagner <kadrlica@fnal.gov>
  More pythonized and streamlined. Added wacky pixel
  finder. Added TODO list...

TODO:
- Find 'dim' columns that aren't quite 'funky'
- Apply other column BPMs before looking for 'wacky' pixels

@author: Alex Drlica-Wagner <kadrlica@fnal.gov>
"""
__author__   = "Alex Drlica-Wagner"
__version__  = "0.0.1"
__revision__ = "$Rev$".strip('$').split()[-1]

import sys
import os
import logging
import copy
import datetime
import subprocess
from collections import OrderedDict as odict

import pyfits
import numpy as np
import fitsio
import scipy.stats
import scipy.ndimage

# Import the maskbits from despyfits:
# imsupport/include/mask_bits.h

# Image geometry
NROW=4096
NCOL=2048
SHAPE=[NROW,NCOL]

# Codes from from $IMSUPPORT_DIR/include/imsupport.h
# Should be able to load with ctypes.CDLL
# I'm ashamed to even think about doing this...
DEFAULT = odict([
        ('BPMDEF_FLAT_MIN' , 1   ), 
        ('BPMDEF_FLAT_MAX' , 2   ), 
        ('BPMDEF_FLAT_MASK', 4   ), 
        ('BPMDEF_BIAS_HOT' , 8   ), 
        ('BPMDEF_BIAS_WARM', 16  ), 
        ('BPMDEF_BIAS_MASK', 32  ), 
        ('BPMDEF_BIAS_COL' , 64  ), 
        ('BPMDEF_EDGE'     , 128 ), 
        ('BPMDEF_CORR'     , 256 ), 
        ('BPMDEF_TAPE_BUMP', 512 ),
        ('BPMDEF_FUNKY_COL', 1024), 
        ('BPMDEF_WACKY_PIX', 2048), 
        ('BPMDEF_GENERIC'  , 4096),  # GENERIC should probably be removed
])

def get_imsupport(basedir=None,filepath='include/imsupport.h'):
    if basedir is None: basedir = os.getenv('IMSUPPORT_DIR') 
    filename = os.path.join(os.getenv('IMSUPPORT_DIR'),filepath)    

    if not os.path.exists(filename): return DEFAULT
        
    ret = odict(DEFAULT)
    for line in open(filename,'r').readlines():
        if not line.startswith('#define BPMDEF_'): continue
        junk,key,value = line.split(' ',3)
        ret[key] = int(value)
    return ret

IMSUPPORT = get_imsupport()

# Just remap for old interface    
BADPIX = odict([
        ('GENERIC'  ,IMSUPPORT['BPMDEF_GENERIC']  ),
        ('EDGE'     ,IMSUPPORT['BPMDEF_EDGE']     ),
        ('CORR'     ,IMSUPPORT['BPMDEF_CORR']     ),
        ('FLAT_MIN' ,IMSUPPORT['BPMDEF_FLAT_MIN'] ),
        ('FLAT_MAX' ,IMSUPPORT['BPMDEF_FLAT_MAX'] ),
        ('BIAS_HOT' ,IMSUPPORT['BPMDEF_BIAS_HOT'] ),
        ('BIAS_WARM',IMSUPPORT['BPMDEF_BIAS_WARM']),
        ('BIAS_COL' ,IMSUPPORT['BPMDEF_BIAS_COL'] ),
        ('FUNKY_COL',IMSUPPORT['BPMDEF_FUNKY_COL']),
        ('WACKY_PIX',IMSUPPORT['BPMDEF_WACKY_PIX']),
        ('TAPE_BUMP',IMSUPPORT['BPMDEF_TAPE_BUMP']),
        ])

BADPIX_DOC = odict([
        ('GENERIC'  , 'Generic (undifferentiated) bad pixel flag'),
        ('EDGE'     , 'Pixels on the edge of the CCD.'   ),
        ('CORR'     , 'Correctable pixels (usually downstream of hot pixels).' ),
        ('FLAT_MIN' , 'Pixels that are dull in the flats.'  ),
        ('FLAT_MAX' , 'Pixels that are hot in the flats.'  ),
        ('BIAS_HOT' , 'Pixels that are hot in the biases.'  ),
        ('BIAS_WARM', 'Pixels that are warm in the biases.'  ),
        ('BIAS_COL' , 'Pixels that are downstream of a hot pixel in the bias.' ),
        ('FUNKY_COL', 'Columns with charge redistribution in sky exposures.' ),
        ('WACKY_PIX', 'Outliers in stacked sky exposures.'  ),
        ('TAPE_BUMP', 'Pixels that reside in tape bumps.' ),
        ])


### Original Marriner codes
#BADPIX_GENERIC  =1
#BADPIX_FLAT_MIN =2
#BADPIX_FLAT_MAX =4
#BADPIX_FLAT_MASK=8
#BADPIX_BIAS_HOT =16
#BADPIX_BIAS_WARM=32
#BADPIX_BIAS_MASK=64
#BADPIX_BIAS_COL =128
#BADPIX_CORR     =256
#BADPIX_EDGE     =512
#BADPIX_FUNKY_COL=1024
#BADPIX_FUNKY_PIX=2048
#BADPIX_TAPE_BUMP=4096

###
### Some general utility functions
###

def positive_int(value):
    ivalue = int(value)
    if ivalue < 0:
        raise argparse.ArgumentTypeError("invalid positive int value: '%s'"%ivalue)
    return ivalue

def idl_histogram(data,bins=None):
    """
    Bins data using numpy.histogram and calculates the
    reverse indices for the entries like IDL.
    
    Parameters:
    data  : data to pass to numpy.histogram
    bins  : bins to pass to numpy.histogram
    Returns:
    hist  : bin content output by numpy.histogram
    edges : edges output from numpy.histogram
    rev   : reverse indices of entries in each bin
 
    Using Reverse Indices:
        h,e,rev = histogram(data, bins=bins)
        for i in range(h.size):
            if rev[i] != rev[i+1]:
                # data points were found in this bin, get their indices
                indices = rev[ rev[i]:rev[i+1] ]
                # do calculations with data[indices] ...
    """
    hist, edges = np.histogram(data, bins=bins)
    digi = np.digitize(data.flat,bins=np.unique(data)).argsort()
    rev = np.hstack( (len(edges), len(edges) + np.cumsum(hist), digi) )
    return hist,edges,rev

def smooth(x,window_len=11,window='hanning'):
    """
    Smoothing function. (From Peter Melchior)
    """
    if x.ndim != 1:
        raise ValueError("Only accepts 1-dimensional arrays.")
    if x.size < window_len:
        raise ValueError("Input vector must be larger than window size.")
    if window_len<3:
        return x
    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError("Window must be: 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")
    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')
    y=np.convolve(w/w.sum(),s,mode='valid')
    return y[(window_len/2):-(window_len/2)]


class AttributeDict(odict):
    # ADW: Careful, memory leak in <=2.7.3
    # http://bugs.python.org/issue1469629    
    def __init__(self, *args, **kwargs):
        super(AttributeDict, self).__init__(*args, **kwargs)
        self.__dict__ = self    
adict = AttributeDict


class BadPixelMasker(object):
    """
    The primary class for creating a bad pixel mask (BPM).
    Defaults are set in the object and passed to the
    argument parser through the 'argparser' method.
    Creation of the BPM occurs via the 'run' method.
    """
    defaults = odict([
            # File IO
            ('outfile',dict(type=str,metavar='bpm.fits',required=True,
                            help='Output bad pixel mask file.'         )), 
            ('flatcor',dict(type=str,metavar='flatcor.lst',
                            help='List of input precal flat files.'    )), 
            ('biascor',dict(type=str,metavar='flatcor.lst',
                            help='List of input precal bias files.'    )), 
            ('images' ,dict(type=str,metavar='flatcor.lst',
                            help='List of input science images.'       )), 

            # Explicitly set pixels
            ('badpix' ,dict(default=None,metavar='badpix.fits',
                            help="Bad pixel map to combine with output BPM"    )),
            ('funkycol',dict(default=None,metavar='funky_pixels.lst',
                             help="Funky column map to combine with output BPM.")),
            ('wackypix',dict(default=None,metavar='wacky_pixels.lst',
                             help="Wacky pixel map to combine with output BPM." )),

            # Tuneable parameters for algorithms
            ('edgesize',dict(default=15,type=positive_int,
                             help="Number of edge pixels to mask."              )),
            ('flatmax',dict(default=1.5,type=float,
                            help="Threshold for hot pixels in the flats."       )),
            ('numfmax',dict(default=1,type=positive_int,
                            help="Number of flats where pixel must be hot."     )),
            ('flatmin',dict(default=0.5,type=float,
                            help="Threshold for cold pixels in the flats."      )),
            ('numfmin',dict(default=1,type=positive_int,
                            help="Number of flats where pixel must be cold."    )),
            ('biashot',dict(default=500.,type=float,
                            help="Threshold for hot pixels in the biases."      )),
            ('numbhot',dict(default=1,type=positive_int,
                            help="Number of biases where pixel must be hot."    )),
            ('biaswarm',dict(default=20.,type=float,
                             help="Threshold for warm pixels in the biases."    )),
            ('numbwarm',dict(default=1,type=positive_int,
                             help="Number of biases where pixel must be warm."  )),
            ('funkymin',dict(default=-5,type=float,
                             help="Threshold for funky column."                 )),
            ('wackynpix',dict(default=4,type=positive_int,                       
                              help="Minimum size of funky pixel group."         )),
            ('wackynsig',dict(default=7,type=float,
                              help="Sigma clip for bad pixels."                 )),

            # Some general properties
            ('generic',dict(action='store_true',                                 
                            help="Set all BPM bits to BADPIX=1"                 )),
            ('badccd',dict(action='store_true',                                  
                           help="Produce trivial BPM for bad CCDs."             )),
            ('outdir',dict(default=None,type=str,
                            help="Output directory for QA."                     )),

        ])

    def __init__(self, options):
        self.opts = self.parse_options(self.defaults, options)
        self.ccdnum = None
        # Currently, outdir must exist
        self.outdir = self.opts.outdir
        if self.outdir and not os.path.exists(self.outdir):
            logging.info("Creating QA directory: %s"%self.outdir)
            os.makedirs(self.outdir)
            #msg = "Debug directory does not exist: %s"%self.outdir
            #raise Exception(msg)
        self.debug = self.outdir is not None

    @staticmethod
    def parse_options(defaults, options):
        opts = adict(defaults)
        for k,v in options.items():
            if k not in defaults:
                msg = "Option '%s' not found in defaults"%k
                raise KeyError(msg)
            opts[k] = v
        return opts

    @classmethod
    def argparser(cls, **kwargs):
        """ 
        Fill ArgumentParser with class defaults as arguments.
        """
        import argparse
        description = "Script for generating bad pixel masks (BPMs) from sets"
        description+= " of flats, biases, and science images."

        kwargs.setdefault('description',description)
        kwargs.setdefault('formatter_class',argparse.ArgumentDefaultsHelpFormatter)
        parser = argparse.ArgumentParser(**kwargs)
        for key,value in cls.defaults.items():
            parser.add_argument('--%s'%key,dest=key,**value)
        return parser

    @staticmethod
    def file_info(fits):
        for i in range(len(fits)):
            try:                   
                header = fits[i].header
            except AttributeError: 
                header = fits[i].read_header()
            if header.get('DES_EXT').strip() == 'IMAGE':
                ext=i
                break
        ccdnum = header.get('CCDNUM')
                
        return ccdnum,ext

    @staticmethod
    def load_file_list(filename):
        """
        Create an array of file names from single file
        or input file list.
     
        Parameters:
          filename : Input FITS file or list of FITS files
     
        Returns:
          filelist : Array of file names
        """
        logging.info("Loading file list: %s"%filename)

        try:
            # Fits file, return as list
            pyfits.open(filename)
            return np.array([filename])
        except :
            # List of files, return list
            return np.loadtxt(filename,dtype=str)
        raise IOError("Failed to load %s"%filename)
    
    @staticmethod
    def load_pixel_list(filename, ccdnum):
        bpm = np.zeros(SHAPE,dtype=int) 
        badpix = np.loadtxt(filename,dtype=int)
        if badpix.shape[-1] < 5:
            raise Exception("Wrong number of columns in %s"%filename)
        if badpix.shape[-1] == 5:
            badpix = np.append(badpix,np.ones(len(badpix)),axis=1)

        for i,xmin,xmax,ymin,ymax,bit in badpix:
            if i != ccdnum: continue
            bpm[ymin-1:ymax,xmin-1:xmax] |= bit
        return bpm

    @staticmethod
    def load_bpm(filename):
        f = pyfits.open(filename)
        bpm = copy.copy(f[0].data)
        return bpm.astype(np.uint16)

    @staticmethod
    def diagnostic(bpm, opts):
        """
        Count bad pixel type and print diagnostic output.
     
        Parameters:
          bpm  : Bad pixel mask array
        Returns:
        """
        logging.info("")
        logging.info("BAD PIXEL MASK SUMMARY")

        badpix =(bpm&~BADPIX['CORR']>0).sum()
        corrpix=(bpm&BADPIX['CORR']>0).sum()
     
        logging.info("Flat low=%i (%i < %g)"%((bpm&BADPIX['FLAT_MIN']>0).sum(),opts.numfmin,opts.flatmin))
        logging.info("Flat high=%i (%i > %g)"%((bpm&BADPIX['FLAT_MAX']>0).sum(),opts.numfmax,opts.flatmax))
        #logging.info("Flat mask=%i"%(bpm&BADPIX_FLAT_MASK>0).sum())
        logging.info("Bias hot=%i (%i > %g)"%((bpm&BADPIX['BIAS_HOT']>0).sum(),opts.numbhot,opts.biashot))
        logging.info("Bias warm=%i (%i > %g)"%((bpm&BADPIX['BIAS_WARM']>0).sum(),opts.numbwarm,opts.biaswarm))
        #logging.info("Bias mask=%i"%(bpm&BADPIX_BIAS_MASK>0).sum())
        logging.info("Bias column=%i"%(bpm&BADPIX['BIAS_COL']>0).sum())
        logging.info("Edge=%i (%i)"%((bpm&BADPIX['EDGE']>0).sum(),opts.edgesize))
        logging.info("Correctable=%i"%(bpm&BADPIX['CORR']>0).sum())
        logging.info("Funky columns=%i"%(bpm&BADPIX['FUNKY_COL']>0).sum())
        logging.info("Wacky pixels=%i"%(bpm&BADPIX['WACKY_PIX']>0).sum())
        logging.info("Generic bad pixels=%i"%(bpm&BADPIX['GENERIC']>0).sum())
     
        badfract  = (100.0*badpix)/bpm.size
        corrfract = (100.0*corrpix)/bpm.size
        logging.info("Total number of uncorrectable pixels masked = %d (%.2f%%)"%(badpix,badfract))
        logging.info("Total number of correctable pixels masked = %d (%.2f%%)"%(corrpix,corrfract))
        logging.info("")

    @staticmethod
    def write_bpm(bpm,filename,ccdnum):
        """
        Write the bad column mask.
     
        Parameters:
          bpm      : Bad pixel mask array
          filename : Output fits file
          ccdnum   : CCD number
     
        Returns:
        """
        dirname = os.path.dirname(filename)
        if dirname and not os.path.exists(dirname):
            msg = "Outfile path does not exist: %s"%filename
            raise IOError(msg)
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
        return hdu

    def write_funky_columns(self, bpm):
        logging.info("Writing funky columns...")
        outfile = os.path.join(self.outdir,"funky_columns_%02d.fits"%self.ccdnum)
        self.write_bpm(bpm,outfile,self.ccdnum)
        subprocess.call('gzip -f %s'%outfile,shell=True)

    def write_wacky_pixels(self, bpm):
        logging.info("Writing wacky pixels...")
        outfile = os.path.join(self.outdir,"wacky_pixels_%02d.fits"%self.ccdnum)
        self.write_bpm(bpm,outfile,self.ccdnum)
        subprocess.call('gzip -f %s'%outfile,shell=True)
        
    def write_median_sky(self, median, header=None):
        logging.info("Writing median sky image...")
        outfile = os.path.join(self.outdir,"median_%02d.fits.fz"%self.ccdnum)
        compress  = 'RICE'
        tile_dims = [1,NCOL]
        ofits = fitsio.FITS(outfile,'rw',clobber=True)
        ofits.write(median,header=header,compress=compress,tile_dims=tile_dims)
        ofits.close()

    def write_funky_devs(self, nites, meds, devs):
        logging.info("Writing column deviations...")
        outfile=os.path.join(self.outdir,"funky_devs_%02d.fits"%self.ccdnum)
        hdu = pyfits.PrimaryHDU(devs)
        col1 = pyfits.Column(name='NITE', format='J', array=nites)
        col2 = pyfits.Column(name='MEDIAN',format='E', array=meds)
        tbhdu = pyfits.new_table(pyfits.ColDefs([col1, col2]))
        hdulist = pyfits.HDUList([hdu, tbhdu])
        hdulist.writeto(outfile,clobber=True)

    def write_funky_list(self,columns):
        logging.info("Writing funky column list...")
        outfile = os.path.join(self.outdir,"funky_%02d.lst"%self.ccdnum)
        out = open(outfile,'w')
        out.write("# CHIP\tCOLUMN\tDEVIATION\n")
        np.savetxt(out,columns,fmt='%-8i%-8i%-8g')
        out.close()


    def run(self,opts=None):
        """
        Create the BPM.
        """
        if opts is None: opts = self.opts
        if opts.badccd: bpm = self.run_bad_ccd(opts)
        else:           bpm = self.run_good_ccd(opts)

        self.diagnostic(bpm,self.opts)
        self.write_bpm(bpm,self.opts.outfile,self.ccdnum)

        return bpm

    def run_bad_ccd(self,opts):
        """
        Create the BPM for a bad CCD.
        """
        flatname = self.load_file_list(opts.flatcor)[0]
        self.ccdnum = self.file_info(pyfits.open(flatname))[0]
        bpm = np.zeros(SHAPE,dtype=int) 
        bpm |= self.mask_edges(opts)
        return bpm

    def run_good_ccd(self,opts):
        """
        Create the BPM for a good CCD.
        """
        bpm = np.zeros(SHAPE,dtype=int) 

        bpm |= self.mask_edges(opts)

        if opts.flatcor:
            bpm |= self.mask_flats(opts)
        if opts.biascor:
            bpm |= self.mask_biases(opts)
        if opts.images:
            bpm |= self.mask_images(opts)
            
        # Clear warm pixels in correctable pixels (bad columns)
        bpm &= ~(BADPIX['BIAS_WARM']*(bpm&BADPIX['CORR']>0))
        # Also clear warm pixesl at edges
        bpm &= ~(BADPIX['BIAS_WARM']*(bpm&BADPIX['EDGE']>0))

        if opts.badpix:
            bpm |= self.mask_badpix(opts)

        # Set generic bad pixel bit
        if opts.generic: 
            bpm = (bpm > 0)*BADPIX['GENERIC']
            
        return bpm

    def mask_badpix(self,opts):
        bpm = np.zeros(SHAPE,dtype=int) 
        logging.info("Loading bad pixels: %s"%opts.badpix)
        if os.path.splitext(opts.badpix)[1]=='.fits':
            bpm |= self.load_bpm(opts.badpix)
        else:
            bpm |= self.load_pixel_list(opts.badpix, self.ccdnum)
        return bpm
        
    def mask_edges(self,opts):
        bpm = np.zeros(SHAPE,dtype=int) 
        logging.info("Masking %i pixel edge"%opts.edgesize)
        if not opts.edgesize: return bpm
        bpm |= BADPIX['EDGE']
        bpm[opts.edgesize:-opts.edgesize,opts.edgesize:-opts.edgesize] &= ~BADPIX['EDGE']
        return bpm

    def mask_flats(self, opts):
        bpm = np.zeros(SHAPE,dtype=int) 
        flatmax_bpm = copy.copy(bpm)
        flatmin_bpm = copy.copy(bpm)

        flatnames = self.load_file_list(opts.flatcor)
        if not len(flatnames): return bpm

        logging.info("Processing %i flats"%len(flatnames))
        for i,flatname in enumerate(flatnames):
            flat = pyfits.open(flatname)
            ccdnum, ext = self.file_info(flat)
            if self.ccdnum is None: self.ccdnum = ccdnum
            if ccdnum != self.ccdnum:
                raise Exception("CCD number does not match")
            image = flat[ext].data
            logging.debug("Processing flat: %s"%flatname)
            flatmax_bpm += (image > opts.flatmax)
            flatmin_bpm += (image < opts.flatmin)

            # Flag masks
            # if opts.flag_mask:
            #     bpm |= flat['MASK'].data

        #Flag bad pixels based on majority logic
        bpm |= ((flatmin_bpm >= opts.numfmin) & (bpm&BADPIX['EDGE']==0))*BADPIX['FLAT_MIN']
        bpm |= (flatmax_bpm >= opts.numfmax)*BADPIX['FLAT_MAX']

        return bpm

    def mask_biases(self,opts):
        bpm = np.zeros(SHAPE,dtype=int) 
        biashot_bpm = copy.copy(bpm)
        biaswarm_bpm = copy.copy(bpm)

        biasnames = self.load_file_list(opts.biascor)
        if not len(biasnames): return bpm

        logging.info("Processing %i biases"%len(biasnames))
        for i,biasname in enumerate(biasnames):
            bias = pyfits.open(biasname)
            ccdnum, ext = self.file_info(bias)
            if self.ccdnum is None: self.ccdnum = ccdnum
            if ccdnum != self.ccdnum:
                raise Exception("CCD number does not match")
         
            image = bias[ext].data
            logging.debug("Processing bias: %s"%biasname)
         
            # Check that the adjacent pixels are not also hot
            nhotbias = (image > opts.biashot)
            nhotbias[1:,:]  *= (image[1:,:]>0.2*image[:-1,:])
            nhotbias[:-1,:] *= (image[:-1,:]>0.2*image[1:,:])
         
            biashot_bpm += nhotbias
            biaswarm_bpm += (image > opts.biaswarm)

        bpm |= (biaswarm_bpm >= opts.numbwarm)*BADPIX['BIAS_WARM']
        bpm |= (biashot_bpm >= opts.numbhot)*BADPIX['BIAS_HOT']

        iy,ix = np.nonzero(bpm & BADPIX['BIAS_HOT'])
        # Flag the surrounding pixels as warm
        for _ix,_iy in zip(ix,iy):
            logging.debug('hotpixel x=%i y=%i',_ix,_iy)
            ymin = _iy-5 if _iy-5 > 0 else 0
            ymax = _iy+6 if _iy+6 < SHAPE[0] else SHAPE[0]
            xmin = _ix-1 if _ix-1 > 0 else 0
            xmax = _ix+2 if _ix+2 < SHAPE[1] else SHAPE[1]
            bpm[ymin:ymax,xmin:xmax] |= BADPIX['BIAS_WARM']
         
            if (ix==_ix).sum() == 1:
                logging.debug("Column=%i is correctable.",_ix)
                # Flag upstream 60 pixel cross talk and downstream
                # excess (including knob?) as correctable
                if self.ccdnum <= 31:
                    ystart = _iy+60
                    if ystart < SHAPE[0]:
                        bpm[ystart::60,_ix] &= ~BADPIX['BIAS_WARM']
                        bpm[ystart::60,_ix] |= BADPIX['CORR']
                    bpm[0:ymin,_ix] &= ~BADPIX['BIAS_WARM']
                    bpm[0:ymin,_ix] |= BADPIX['CORR']
                else:
                    ystart = _iy-60
                    if ystart >= 0:
                        bpm[ystart::-60,_ix] &= ~BADPIX['BIAS_WARM']
                        bpm[ystart::-60,_ix] |= BADPIX['CORR']
                    bpm[ymax:,_ix] &= ~BADPIX['BIAS_WARM']
                    bpm[ymax:,_ix] |= BADPIX['CORR']
            else:
                logging.debug("Column=%i is NOT correctable.",_ix)
                bpm[:,_ix] |= BADPIX['BIAS_COL']
        return bpm

    def mask_images(self,opts):
        """
        Find funky columns and funky pixels from science images.
        """
        bpm = np.zeros(SHAPE,dtype=int) 

        imagenames = self.load_file_list(opts.images)
        if not len(imagenames): return bpm

        data = np.zeros([len(imagenames)]+SHAPE, dtype=float)
        if (not opts.funkycol) or (not opts.wackypix):
            logging.info("Processing %i images"%len(imagenames))
            for i,imagename in enumerate(imagenames):
                logging.debug("Processing image: %s"%imagename)
                fits = fitsio.FITS(imagename)
                ccdnum, ext = self.file_info(fits)
                if self.ccdnum is None: self.ccdnum = ccdnum
                header = fits[ext].read_header()
                if header['CCDNUM'] != self.ccdnum:
                    raise Exception("CCD number does not match")
                image = fits[ext].read()
                data[i] = image

        if not opts.funkycol:
            columns = self.find_funky_columns(data, opts)
        else:
            logging.info("Loading funky columns: %s"%opts.funkycol)
            columns = np.loadtxt(opts.funkycol)
        funky_bpm = self.find_badcol(columns,self.ccdnum,opts.funkymin)
        if self.debug: self.write_funky_columns(funky_bpm)

        if not opts.wackypix:
            logging.info("Calculating median sky image...")
            median = np.median(data,axis=0)

            if self.debug: 
                header = fitsio.read_header(imagenames[0],ext)
                self.write_median_sky(median,header)

            wacky_bpm = self.find_wacky_pixels(median, opts)
        else:
            logging.info("Loading wacky pixels: %s"%opts.wackypix)
            wacky_bpm = self.load_pixel_list(opts.wackypix, self.ccdnum)

        wacky_bpm *= (funky_bpm == 0)
        logging.info("Number of wacky pixels: %s"%(wacky_bpm>0).sum())
        if self.debug: self.write_wacky_pixels(wacky_bpm)
                        
        bpm |= funky_bpm*BADPIX['FUNKY_COL']
        bpm |= wacky_bpm*BADPIX['WACKY_PIX']
        return bpm

    def find_funky_columns(self, data, opts):
        """
        Find funky columns.
        """
        logging.info("Calculating column deviations...")
        EDGE = opts.edgesize

        nites = np.arange(data.shape[0])
        devs = np.zeros((data.shape[0],NCOL), dtype=float)
        meds = np.zeros(data.shape[0], dtype=float)

        devs[:,EDGE:-EDGE] = np.median(data,axis=1)[:,EDGE:-EDGE]
        for i,d in enumerate(data):
            meds[i] = np.median(d)
            devs[i][EDGE:-EDGE] -= meds[i]
            #devs[i][EDGE:-EDGE] = np.median(d,axis=0)[EDGE:-EDGE] - meds[i]
         
        if self.debug: self.write_funky_devs(nites,meds,devs)
        
        logging.info("Finding funky columns...")
        stack = np.median(devs,axis=0)
        filtered = stack.copy()
        filtered[np.fabs(stack) > 10] = 0
        smoothed = np.concatenate((smooth(filtered[EDGE:NCOL/2],101),
                                   smooth(filtered[NCOL/2:-EDGE],101)))
        stack[EDGE:-EDGE] -= smoothed
        funky_cols = np.nonzero(abs(stack)>1)[0]
        columns = np.array([self.ccdnum*np.ones(len(funky_cols)),funky_cols,stack[funky_cols]]).T

        if self.debug: self.write_funky_list(columns)

        bpm = self.find_badcol(columns,self.ccdnum,opts.funkymin)
        return bpm

    @staticmethod
    def find_badcol(columns,ccdnum,threshold=-5):
        """
        Mask 'funky' bad columns.
        NOTE: This was translated from C and could probably be pythonized.

        Parameters:
        columns     : Input bad column list
        ccdnum   : CCD being examined
        Returns:
        mask     : Bad column mask array
        """
     
        bpm = np.zeros(SHAPE,dtype=np.uint16)
     
        #columns = np.loadtxt(filename)
        # Empty file
        if columns.size == 0: 
            logging.info("No bad columns found.")
            return bpm
        elif columns.ndim == 1: 
            columns = np.array([columns])
     
        columns = columns[columns[:,0] == ccdnum]
        dev = np.zeros(NCOL)
        bcol = np.zeros(NCOL)
        seq = 0; tot = 0
     
        for ccd,col,d in columns:
            logging.debug("ccd=%i col=%i dev=%f"%(ccd,col,d))
            dev[col] = d
        
        for i in np.arange(NCOL/2):
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
                logging.debug("start=%i size=%f sum=%f"%(start,size,tot))
     
            #See if start of a new sequence
            seq = 0
            if (dev[i]>threshold): continue
            seq = i
            start = seq
            size = dev[i]
            tot = size
            bcol[i] = True
     
        if (seq!=0 and np.fabs(tot)/size<1.0):
            logging.debug("start=%i size=%f sum=%f"%(start,size,tot))
     
        seq = 0
        for i in np.arange(NCOL/2,NCOL)[::-1]:
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
                logging.debug("start=%i size=%f sum=%f"%(start,size,tot))
     
            #See if start of a new sequence
            seq = 0;
            if (dev[i]>threshold): continue;
            seq = i
            start = seq
            size = dev[i]
            tot = size
            bcol[i] = True
     
        if (seq!=0 and np.fabs(tot)/size<1.0):
            logging.debug("start=%i size=%f sum=%f"%(start,size,tot))
     
        logging.info("Found %i bad columns"%bcol.sum());
        for i in np.nonzero(bcol)[0]:
            logging.debug("Masking bad column=%i"%i)        
        bpm[:,np.nonzero(bcol)[0]] |= True
     
        return bpm

    def find_wacky_pixels(self,median,opts):
        """
        Find islands of wacky pixels.
        """
        logging.info("Finding wacky pixels...")
        EDGE = opts.edgesize
        clipped = median[EDGE:-EDGE,EDGE:-EDGE]
        c,cmin,cmax = scipy.stats.sigmaclip(clipped,opts.wackynsig,opts.wackynsig)
        logging.info("Median sigma-clipped at: %s %s"%(cmin,cmax))
         
        searchim = np.zeros(SHAPE,dtype=int)
        searchim[EDGE:-EDGE,EDGE:-EDGE] = (clipped<cmin)|(clipped>cmax)
        structure = scipy.ndimage.generate_binary_structure(2,2)
        label,nlabel = scipy.ndimage.label(searchim,structure)
        hist,edges,rev = idl_histogram(label,bins=np.arange(nlabel+2))
        good, = np.where((hist >= opts.wackynpix) & (edges[:-1] != 0))
         
        bpm = np.zeros(SHAPE,dtype=np.uint16)
        nislands = 0
        for i in good:
            i1a=rev[rev[i]:rev[i+1]]
            xpix = i1a % NCOL
            ypix = i1a / NCOL
            if xpix.ptp() == 0:
                logging.debug("Ignoring island at: %s,%s"%(xpix.mean(),ypix.mean()))
                continue
            if ypix.ptp() > NCOL//8:
                logging.debug("Ignoring island at: %s,%s"%(xpix.mean(),ypix.mean()))
                continue
            if ypix.ptp() > NCOL//16 and xpix.ptp() < 2:
                logging.debug("Ignoring island at: %s,%s"%(xpix.mean(),ypix.mean()))
                continue
            bpm[ypix,xpix] = True
            nislands += 1
         
        # Don't grab artifacts of bad columns
        logging.info("Number of wacky islands: %s"%nislands)
        return bpm

if __name__ == "__main__":
    parser = BadPixelMasker.argparser()
    parser.add_argument('-v','--verbose',default=2,choices=range(4),type=int,
                        help='Verbosity level')
    opts = parser.parse_args()

    verbose = vars(opts).pop('verbose')
    logging.basicConfig(level=10*(4-verbose),
                        format='MKBPM:%(levelname)s: %(message)s')

    masker = BadPixelMasker(vars(opts))
    bpm    = masker.run()
