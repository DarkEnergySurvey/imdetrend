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
import fitsio
import subprocess

import scipy.stats
import scipy.ndimage

# Image geometry
NROW=4096
NCOL=2048
SHAPE=[NROW,NCOL]

#generic bad pixel code
#BADPIX=1
#test specific bad pixel codes
#ADW: BADPIX should be different from BADPIX_FLAT_MIN
BADPIX_FLAT_MIN=1
BADPIX_FLAT_MAX=2
BADPIX_FLAT_MASK=4
BADPIX_BIAS_HOT=8
BADPIX_BIAS_WARM=16
BADPIX_BIAS_MASK=32
BADPIX_BIAS_COL=64
BADPIX_EDGE=128
BADPIX_CORR=256
BADPIX_FUNKY_COL=512
BADPIX_FUNKY_PIX=1024
BADPIX=2048

class BadPixelMasker(object):
    #ADW: Should break this up a bit more into
    #file loading functions and bpm processing
    #functions
    def __init__(self,*args,**kwargs):
        self.ccdnum = None
        self.outdir = opts.outdir # outdir must exist
        #if opts.debug and not os.path.exists(self.outdir):
        #    os.makedirs(self.outdir)

    @staticmethod
    def load(filename):
        """
        Create an array of file names from single file
        or input file list.
     
        Parameters:
          filename : Input FITS file or list of FITS files
     
        Returns:
          filelist : Array of file names
        """
        try:
            # Fits file, return as list
            pyfits.open(filename)
            return np.array([filename])
        except :
            # List of files, return list
            return np.loadtxt(filename,dtype=str)
        raise IOError("Failed to load %s"%filename)

    @staticmethod
    def load_bpm(filename):
        f = pyfits.open(filename)
        bpm = copy.copy(f[0].data)
        return bpm.astype(np.uint16)

    @staticmethod
    def info(fits):
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


    @staticmethod
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


    @staticmethod
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

    @staticmethod
    def badcol(data,ccdnum,threshold=-5):
        """
        Mask 'funky' bad columns.
        NOTE: This was translated from C and could probably be pythonized.

        Parameters:
        data     : Input bad column list
        ccdnum   : CCD being examined
        Returns:
        mask     : Bad column mask array
        """
     
        bpm = np.zeros(SHAPE,dtype=np.uint16)
     
        #data = np.loadtxt(filename)
        # Empty file
        if data.size == 0: 
            logging.info("No bad columns found.")
            return bpm
        elif data.ndim == 1: 
            data = np.array([data])
     
        data = data[data[:,0] == ccdnum]
        dev = np.zeros(NCOL)
        bcol = np.zeros(NCOL)
        seq = 0; tot = 0
     
        for ccd,col,d in data:
            logging.debug("ccd=%i col=%i dev=%f"%(ccd,col,d))
            dev[col] = d
        
        for i in np.arange(NCOL//2):
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

    def run_badccd(self,opts):
        flatname = self.load(opts.flatcor)[0]
        self.ccdnum = self.info(pyfits.open(flatname))[0]
        bpm = np.zeros(SHAPE,dtype=int) 
        bpm |= self.mask_edges(opts)
        return bpm

    def run_all(self,opts):
        bpm = np.zeros(SHAPE,dtype=int) 
        bpm |= self.run_flats(opts)
        bpm |= self.run_biases(opts)
        bpm |= self.run_images(opts)
        
        if opts.edgesize: 
            bpm |= self.mask_edges(opts)
            #logging.info("Masking %i pixel edge"%opts.edgesize)
            #bpm |= BADPIX_EDGE
            #bpm[opts.edgesize:-opts.edgesize,opts.edgesize:-opts.edgesize] &= ~BADPIX_EDGE

        # Clear warm bias pixels in correctable pixels (bad columns)
        bpm &= ~(BADPIX_BIAS_WARM*(bpm&BADPIX_CORR>0))
        # Also at edges
        bpm &= ~(BADPIX_BIAS_WARM*(bpm&BADPIX_EDGE>0))

        if opts.badpix:
            bpm |= self.mask_badpix(opts)

        # Set generic bad pixel bit
        if opts.generic:
            #bpm = (bpm > 0)*BADPIX
            bpm = (bpm > 0)
        return bpm

    def run_flats(self, opts):
        bpm = np.zeros(SHAPE,dtype=int) 
        flatmax_bpm = copy.copy(bpm)
        flatmin_bpm = copy.copy(bpm)

        flatnames = self.load(opts.flatcor)
        logging.info("Using flatcor (value>%.2f and value<%.2f) and biascor (>%.2f) to identify bad pixels",opts.flatmax,opts.flatmin,opts.biashot)

        logging.info("Processing %i flats"%len(flatnames))
        for i,flatname in enumerate(flatnames):
            flat = pyfits.open(flatname)
            ccdnum, ext = self.info(flat)
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
        bpm |= ((flatmin_bpm >= opts.numfmin) & (bpm&BADPIX_EDGE==0))*BADPIX_FLAT_MIN
        bpm |= (flatmax_bpm >= opts.numfmax)*BADPIX_FLAT_MAX

        return bpm

    def run_biases(self,opts):
        bpm = np.zeros(SHAPE,dtype=int) 
        biashot_bpm = copy.copy(bpm)
        biaswarm_bpm = copy.copy(bpm)

        biasnames = self.load(opts.biascor)
        
        logging.info("Processing %i biases"%len(biasnames))
        for i,biasname in enumerate(biasnames):
            bias = pyfits.open(biasname)
            ccdnum, ext = self.info(bias)
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

        bpm |= (biaswarm_bpm >= opts.numbwarm)*BADPIX_BIAS_WARM
        bpm |= (biashot_bpm >= opts.numbhot)*BADPIX_BIAS_HOT

        iy,ix = np.nonzero(bpm & BADPIX_BIAS_HOT)
        # Flag the surrounding pixels as warm
        for _ix,_iy in zip(ix,iy):
            logging.debug('hotpixel x=%i y=%i',_ix,_iy)
            ymin = _iy-5 if _iy-5 > 0 else 0
            ymax = _iy+6 if _iy+6 < SHAPE[0] else SHAPE[0]
            xmin = _ix-1 if _ix-1 > 0 else 0
            xmax = _ix+2 if _ix+2 < SHAPE[1] else SHAPE[1]
            bpm[ymin:ymax,xmin:xmax] |= BADPIX_BIAS_WARM
         
            if (ix==_ix).sum() == 1:
                logging.debug("Column=%i is correctable.",_ix)
                # Flag upstream 60 pixel cross talk and downstream
                # excess (including knob?) as correctable
                if self.ccdnum <= 31:
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
                    bpm[ymax:,_ix] &= ~BADPIX_BIAS_WARM
                    bpm[ymax:,_ix] |= BADPIX_CORR
            else:
                logging.debug("Column=%i is NOT correctable.",_ix)
                bpm[:,_ix] |= BADPIX_BIAS_COL
        return bpm

    def run_images(self,opts):
        """
        Find funky columns and funky pixels from science images.
        """
        bpm = np.zeros(SHAPE,dtype=int) 

        EDGE = opts.edgesize
        imagenames = self.load(opts.images)

        data = np.zeros([len(imagenames)]+SHAPE, dtype=float)

        if (not opts.funkycol) or (not opts.funkypix):
            logging.info("Processing %i images"%len(imagenames))
            for i,imagename in enumerate(imagenames):
                logging.debug("Processing image: %s"%imagename)
                fits = fitsio.FITS(imagename)
                ccdnum, ext = self.info(fits)
                if self.ccdnum is None: self.ccdnum = ccdnum
                header = fits[ext].read_header()
                if header['CCDNUM'] != self.ccdnum:
                    raise Exception("CCD number does not match")
                image = fits[ext].read()
                data[i] = image

        if not opts.funkycol:
            column_bpm = self.mask_funky_columns(data, opts)
        else:
            logging.info("Loading funky columns: %s"%opts.funkycol)
            output = np.loadtxt(opts.funkycol)
            column_bpm = self.badcol(output,self.ccdnum,opts.funkymin)

        if opts.debug:
            logging.info("Writing funky columns...")
            outfile = os.path.join(self.outdir,"funky_column_%02d.fits"%self.ccdnum)
            ofits = fitsio.FITS(outfile,'rw',clobber=True)
            ofits.write(column_bpm)
            ofits.close()
            subprocess.call('gzip -f %s'%outfile,shell=True)

        if not opts.funkypix:
            logging.info("Calculating median sky image...")
            median = np.median(data,axis=0)
             
            if opts.debug:
                logging.info("Writing median sky image...")
                outfile = os.path.join(self.outdir,"median_%02d.fits.fz"%self.ccdnum)
                compress  = 'RICE'
                tile_dims = [1,NCOL]
                ifits = fitsio.FITS(imagenames[0])
                ofits = fitsio.FITS(outfile,'rw',clobber=True)
                ofits.write(median,header=ifits[ext].read_header(),compress=compress,tile_dims=tile_dims)
                ofits.close()

            pixel_bpm = self.mask_funky_pixels(median, opts)
        else:
            logging.info("Loading funky pixels: %s"%opts.funkypix)
            output = np.loadtxt(opts.funkypix)
            pixel_bpm = np.zeros(SHAPE,dtype=np.uint16)
            for i,xmin,xmax,ymin,ymax in output:
                if i != self.ccdnum: continue
                pixel_bpm[ymin-1:ymax,xmin-1:xmax] = True

        pixel_bpm *= (column_bpm == 0)
        logging.info("Number of funky pixels: %s"%pixel_bpm.sum())

        if opts.debug:
            logging.info("Writing funky pixels...")
            outfile = os.path.join(self.outdir,"funky_pixels_%02d.fits"%self.ccdnum)
            ofits = fitsio.FITS(outfile,'rw',clobber=True)
            ofits.write(pixel_bpm)
            ofits.close()
            subprocess.call('gzip -f %s'%outfile,shell=True)
                        
        bpm |= column_bpm*BADPIX_FUNKY_COL
        bpm |= pixel_bpm*BADPIX_FUNKY_PIX

        return bpm

    # This would be nice...
    def mask_edges(self,data,mask,opts):
        pass
    def mask_funky_columns(self,data,mask,opts):
        pass
    def mask_funky_pixels(self,data,mask,opts):
        pass
    def mask_flats(self,data,mask,opts):
        pass
    def mask_biases(self,data,mask,opts):
        pass

    def mask_badpix(self,opts):
        logging.info("Loading bad pixels: %s"%opts.badpix)
        if os.path.splitext(opts.badpix)[1] == '.fits':
            badpix_bpm = pyfits.open(opts.badpix)[0].data
        else:
            badpix_bpm = np.zeros(SHAPE,dtype=int) 
            badpix = np.loadtxt(opts.badpix)
            for i,xmin,xmax,ymin,ymax in badpix:
                if i != self.ccdnum: continue
            badpix_bpm[ymin-1:ymax,xmin-1:xmax] = True
        return badpix_bpm
        
    def mask_edges(self,opts):
        logging.info("Masking %i pixel edge"%opts.edgesize)
        edge_bpm = np.zeros(SHAPE,dtype=int) 
        edge_bpm |= BADPIX_EDGE
        edge_bpm[opts.edgesize:-opts.edgesize,opts.edgesize:-opts.edgesize] &= ~BADPIX_EDGE
        return edge_bpm

    def mask_funky_columns(self, data, opts):
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
         
        if opts.debug:
            logging.info("Writing column deviations...")
            hdu = pyfits.PrimaryHDU(devs)
            col1 = pyfits.Column(name='NITE', format='J', array=nites)
            col2 = pyfits.Column(name='MEDIAN',format='E', array=meds)
            tbhdu = pyfits.new_table(pyfits.ColDefs([col1, col2]))
            hdulist = pyfits.HDUList([hdu, tbhdu])
            outfile=os.path.join(self.outdir,"funky_devs_%02d.fits"%self.ccdnum)
            hdulist.writeto(outfile,clobber=True)
        
        logging.info("Finding funky columns...")
        stack = np.median(devs,axis=0)
        filtered = stack.copy()
        filtered[np.fabs(stack) > 10] = 0
        smoothed = np.concatenate((self.smooth(filtered[EDGE:NCOL/2],101),
                                   self.smooth(filtered[NCOL/2:-EDGE],101)))
        stack[EDGE:-EDGE] -= smoothed
        funky_cols = np.nonzero(abs(stack)>1)[0]
        output = np.array([self.ccdnum*np.ones(len(funky_cols)),funky_cols,stack[funky_cols]]).T
        column_bpm = self.badcol(output,self.ccdnum,opts.funkymin)        

        if opts.debug:
            logging.info("Writing funky column list...")
            outfile = os.path.join(self.outdir,"funky_%02d.lst"%self.ccdnum)
            out = open(outfile,'w')
            out.write("# CHIP\tCOLUMN\tDEVIATION\n")
            np.savetxt(out,output,fmt='%-8i%-8i%-8g')
            out.close()

        return column_bpm

    def mask_funky_pixels(self,median,opts):
        """
        Find islands of funky pixels.
        """
        logging.info("Finding funky pixels...")
        EDGE = opts.edgesize
        clipped = median[EDGE:-EDGE,EDGE:-EDGE]
        c,cmin,cmax = scipy.stats.sigmaclip(clipped,opts.funkynsig,opts.funkynsig)
        logging.info("Sigma-clipped at: %s %s"%(cmin,cmax))
         
        searchim = np.zeros(SHAPE,dtype=int)
        searchim[EDGE:-EDGE,EDGE:-EDGE] = (clipped<cmin)|(clipped>cmax)
        structure = scipy.ndimage.generate_binary_structure(2,2)
        label,nlabel = scipy.ndimage.label(searchim,structure)
        hist,edges,rev = self.idl_histogram(label,bins=np.arange(nlabel+2))
        good, = np.where((hist >= opts.funkynpix) & (edges[:-1] != 0))
         
        pixel_bpm = np.zeros(SHAPE,dtype=np.uint16)
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
            pixel_bpm[ypix,xpix] = True
            nislands += 1
         
        # Don't grab artifacts of bad columns
        logging.info("Number of funky islands: %s"%nislands)
        return pixel_bpm


def diagnostics(bpm):
    """
    Count bad pixel type and print diagnostic output.

    Parameters:
      bpm  : Bad pixel mask array
    Returns:
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
    logging.info("Funky columns=%i"%(bpm&BADPIX_FUNKY_COL>0).sum())
    logging.info("Funky pixels=%i"%(bpm&BADPIX_FUNKY_PIX>0).sum())

    badfract = (100.0*badpix)/bpm.size
    corrfract = (100.0*corrpix)/bpm.size
    logging.info("Total number of uncorrectable pixels masked = %d (%.2f%%)"%(badpix,badfract))
    logging.info("Total number of correctable pixels masked = %d (%.2f%%)"%(corrpix,corrfract))

def positive_int(value):
    ivalue = int(value)
    if ivalue < 0:
        raise argparse.ArgumentTypeError("invalid positive int value: '%s'"%ivalue)
    return ivalue

if __name__ == "__main__":
    import argparse
    description = "Script for generating bad pixel masks from a set of flats, biases, and science images."
    parser = argparse.ArgumentParser(description=description,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('flatcor',type=str,help='List of precal flat files.')
    parser.add_argument('biascor',type=str,help='List of precal bias files.')
    parser.add_argument('images',type=str,help='List of science images.')
    parser.add_argument('outfile',type=str,help='Output BPM name')

    parser.add_argument('--funkycol',default=None,metavar='funky_columns.lst',
                        help="Funky column map to combine with output BPM.")
    parser.add_argument('--funkypix',default=None,metavar='funky_pixels.lst',
                        help="Funky pixel map to combine with output BPM.")
    parser.add_argument('--badpix',default=None,metavar='bpm.fits',
                        help="Bad pixel map to combine with output BPM")

    parser.add_argument('--flatmax',default=1.5,type=float,
                        help="Threshold for hot pixels in the flats.")
    parser.add_argument('--numfmax',default=1,type=positive_int,
                        help="Number of flats where pixel must be hot.")
    parser.add_argument('--flatmin',default=0.5,type=float,
                        help="Threshold for cold pixels in the flats.")
    parser.add_argument('--numfmin',default=1,type=positive_int,
                        help="Number of flats where pixel must be cold.")
    parser.add_argument('--biashot',default=500.,type=float,
                        help="Threshold for hot pixels in the biases.")
    parser.add_argument('--numbhot',default=1,type=positive_int,
                        help="Number of biases where pixel must be hot.")
    parser.add_argument('--biaswarm',default=20.,type=float,
                        help="Threshold for warm pixels in the biases.")
    parser.add_argument('--numbwarm',default=1,type=positive_int,
                        help="Number of biases where pixel must be warm.")
    parser.add_argument('--funkymin',default=-5,type=float,
                        help="Threshold for funky column.")
    parser.add_argument('--funkynsig',default=7,type=float,
                        help="Sigma clip for bad pixels.")
    parser.add_argument('--funkynpix',default=4,type=positive_int,
                        help="Minimum size of funky pixel group.")
    parser.add_argument('--edgesize',default=15,type=positive_int,
                        help="Number of edge pixels to mask.")
    parser.add_argument('--generic',action='store_true',
                        help="Set all BPM bits to BADPIX=1")
    parser.add_argument('--badccd',action='store_true',
                        help="Produce trivial BPM for bad CCDs.")

    parser.add_argument('-v','--verbose',default=2,choices=range(4),type=int,
                        help='Verbosity level')
    parser.add_argument('--debug',action='store_true',
                        help="Write debug output files")
    parser.add_argument('--outdir',default='debug',
                        help="Path to debug output files")

    opts = parser.parse_args()
    if opts.debug: opts.verbose = 4
    logging.basicConfig(level=10*(4-opts.verbose),
                        format='MKBPM:%(levelname)s: %(message)s')
    
    masker = BadPixelMasker()
    if opts.badccd:
        bpm = masker.run_badccd(opts)
    else:
        bpm = masker.run_all(opts)
    diagnostics(bpm)
    masker.write(bpm,opts.outfile,masker.ccdnum)
