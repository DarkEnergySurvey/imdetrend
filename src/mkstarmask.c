/*$Id$ */
/*subroutine to mask out bright stars using USNOB catalog*/

#include "imsupport.h"
#include <math.h>
static const char *svn_id = "$Id$";

void mkstarmask(input,output,filter,astrostdsfile,horiztrails,exposure,flag_verbose,flag_nointerpolate)
     desimage *input,*output;
     int flag_verbose,horiztrails,flag_nointerpolate;
     char *filter, *astrostdsfile ;
     float exposure;
{
  int	i,j,x,y,k,l,flag_fringe=NO,flag_illum=YES,match,radbins,xpos,ypos,
    minsize=4,maxsize=10,loc,dx,dy,xmin,xmax,ymin,ymax,nmask,
    flag_median=NO,flag_average=NO,count,scalenum,
    nvec,xlen,ylen,totpix,mkpath(),ccdnum=0,ncompare=0,
    flag_illumcor_compare=NO,flag_fringecor_compare=NO,
    scaleregionn[4],ctr,found,*imageval,
    keysexist,numcrays,interp,xcent,ycent,
    pixrow,pixcol,*circlepix ;
  char    event[10000];
	        
  static int status=0;
  long	ranseed=0,seed = -15 ;
  float	ran1(),mode,scalefactor,*scalesort,mean=0,th,logradius,
    goodpixx[10],goodpixy[10],**pixval,outy,erry,radius,mag,
    badrow,badcolumn,ramin,ramax,decmin,decmax,*b2,*r2,*i2,xmask,ymask,krad,
    sigma,*centerx,*centery,*radvalues,*median,gasdev(),*sigimage,
    getmagnitude(),skymean;
  double	value,*tra,*tdec,wcspos[2],rawpos[2],wcsedges[4][2],*circlemean ;
  void	shell(),reportevt(),retrievescale(),extractstars_usnob(),wcsCorners(),pixelhisto(),growcircleradius();
  void    polin2(),creategrid(),bleedingtrails(),readstdstable(),bleedtrailsx(),wcs2xy();	
  time_t	tm;
  float	interp_noise=1.0,interp_fwhm=0.0;
  int	num = 0,goodpix,badpix,n,flag_growcirclesize=1;


  imageval=(int *)calloc(output->npixels,sizeof(int));
  for (n=0;n<output->npixels;n++)imageval[n] =-100;
  if (imageval==NULL) reportevt(flag_verbose,STATUS,5,
				"Calloc for imageval failed");

  /* ************************************ */
  /*   get RA and DEC of image corners 	*/
  /* ************************************ */
  wcsCorners(input->name,wcsedges,0); 
  ramin = wcsedges[1][0]; ramax = wcsedges[3][0] ;
  decmin = wcsedges[1][1]; decmax = wcsedges[3][1] ;
  /* ************************************ */
  /* GET STARS FROM USNOB CATALOG  */
  /* ************************************ */
  tra=(double *)calloc(10000,sizeof(double));
  tdec=(double *)calloc(10000,sizeof(double));
  b2=(float *)calloc(10000,sizeof(float));
  r2=(float *)calloc(10000,sizeof(float));
  i2=(float *)calloc(10000,sizeof(float));
  centerx = (float *)calloc(10000,sizeof(float));
  centery = (float *)calloc(10000,sizeof(float));
  radvalues = (float *)calloc(10000,sizeof(float));
  median = (float *)calloc(10000,sizeof(float));
  sigimage = (float *)calloc(10000,sizeof(float));
  circlemean = (double *)calloc(10000,sizeof(float));
  circlepix = (int *)calloc(10000,sizeof(float));
	
  match = 0 ;

  /*	extractstars_usnob(ramin,ramax,decmin,decmax,tra,tdec,b2,r2,i2,&match) ; */
  readstdstable(astrostdsfile,ramin,ramax,decmin,decmax,tra,tdec,b2,r2,i2,&match);
  sprintf(event,"Returned %d stars from USNOB catalog to be used for masking\n",match);
  reportevt(flag_verbose,STATUS,1,event);
  if (match < 1) {
    sprintf(event,"No bright stars found in catalog");
    reportevt(flag_verbose,STATUS,3,event);
  }
  /* copy the correct output name into place */
  /********************/
  /* APPLY STAR MASK  */
  /********************/
  if (match > 0)
    for (j=0; j<match;j++) 
      {

	wcspos[0] = tra[j],
	  wcspos[1] = tdec[j] ;	
	/* get radius of circle */
	mag = getmagnitude(filter,b2[j],r2[j],i2[j]) ;
		
	/* DC3 only		logradius = -0.3316*mag + 5.2977 ; */
	/* initialize radius to 0 */
	radvalues[j]= 0;
	if (mag>13) continue;
	/* if radius < 100 scale radius with exposure */
	printf("Masking star no %d\n",j);	
	logradius = -0.02382*mag + 2.14406 ; 
	radius = pow(10,logradius) ;
	/*convert RA,DEC to x and y position*/
	wcs2xy(input->name,wcspos,rawpos,flag_verbose); 
	radvalues[j] = radius ; 
	centerx[j] = rawpos[0] ;
	centery[j] = rawpos[1] ;				     	
	radbins = (int) radius ;
	sprintf(event,"x and y positions of circle and radius and magnitude = %f %f %f %f\n",centerx[j],centery[j],radius,mag); 
	reportevt(flag_verbose,STATUS,1,event);		
	for (krad = 0;krad<radius;krad+=0.1)
	  for(th= 0;th<=360;th+=0.25)
	    {
	      xmask = krad*cos(th*3.1416/180.0);
	      ymask = krad*sin(th*3.1416/180.0); 	  		
	      xpos = (int) (xmask + rawpos[0]) ;
	      ypos = (int) (ymask + rawpos[1]) ;	
	      if(xpos<0)xpos=0.0 ;
	      if(ypos<0)ypos=0.0 ;
	      if(xpos >= output->axes[0])xpos=output->axes[0] -1 ;
	      if(ypos >= output->axes[1])ypos=output->axes[1] - 1;
	      i = ypos*output->axes[0] + xpos ;		 
	      if (!(output->mask[i] & BADPIX_STAR))
		{  
		  if (i>output->npixels)
		    {
		      sprintf(event,"Pixel indices out of range for xpos = %d,ypos = %d",xpos,ypos);
		      reportevt(flag_verbose,STATUS,4,event);
		      i = output->npixels;

		    }
		  output->mask[i] |= BADPIX_STAR ;
		  output ->varim[i] =  0.0 ; 
		  nmask++ ;
		  circlemean[j] += output->image[i] ;			  				
		  circlepix[j] += 1;	

		  /*save the USNOB star number in imageval*/
		  imageval[i] = j ;
		}	
	    }
	circlemean[j] = circlemean[j]/circlepix[j] ;

      }

  /*get skymean*/
  ctr = 0;
  skymean=0 ;
  for(num=0;num<output->npixels;num++)
    if (!output->mask[num])
      {skymean+=output->image[num];ctr+=1;}
  skymean /= ctr ;
  /*Find the median values in a section twice that of the masked circle*/
  if (match >0)
    for (j=0;j<match;j++) 
      {
	/*  printf ("before circlemean=%d %f\n",j,radvalues[j]); */ 
	if ((flag_growcirclesize)&&(exposure>100))growcircleradius(output,&radvalues[j],centerx[j],centery[j],skymean,circlemean[j],imageval,j);
	/*    printf ("after circlemean=%d %f\n",j,radvalues[j]); */
	radius = radvalues[j] ;	
	ctr = 0 ;

	scaleregionn[0] = (int) (centerx[j] - 3.0*radius) ;
	scaleregionn[1] = (int) (centerx[j] + 3.0*radius) ;
	scaleregionn[2] = (int) (centery[j] - 3.0*radius) ;
	scaleregionn[3] = (int) (centery[j] + 3.0*radius) ;
	if (scaleregionn[0]<0)scaleregionn[0]= 0;
	if (scaleregionn[1]>output->axes[0])scaleregionn[1]= output->axes[0];
	if (scaleregionn[2]<0)scaleregionn[2]= 0;
	if (scaleregionn[3]>output->axes[1])scaleregionn[3]= output->axes[1];


	scalenum=(scaleregionn[1]-scaleregionn[0]+1)*
	  (scaleregionn[3]-scaleregionn[2]+1);
	if (scalenum > 0)
	  {
	    scalesort=(float *)calloc(scalenum,sizeof(float));
                        
	    if (scalesort==NULL) {
	      reportevt(flag_verbose,STATUS,4,"Calloc of scalesort failed");
	      exit(0); 
	    }
	    retrievescale(output,scaleregionn,scalesort,flag_verbose,&scalefactor,
			  &mode,&sigma) ; 

	    if(fabs(sigma)<7.0)
	      { sprintf (event,"retreivescale sigma too small calling pixelhisto\n");
		reportevt(flag_verbose,STATUS,2,event);
		pixelhisto(output,scaleregionn,scalesort,flag_verbose,&scalefactor,
			   &mode,&sigma);
		if (sigma<7.0)sigma=7.0;} 	      
	    median[j] = scalefactor ; 
	    sigimage[j] = sigma ;

	    free(scalesort);
	    if (circlemean[j] < scalefactor + sigma)circlemean[j] = 0; 
	    sprintf(event,"Circle median and std. deviation= %f %f\n",median[j],sigimage[j]);
	    reportevt(flag_verbose,STATUS,1,event); 
	  }		
      }
  /* fill the bad pixel values with the correct median value corresponding to  */
  /* a rectangular region twice as big as the star */
  for (i=0;i<output->npixels;i++) 
    { 
      if (output->mask[i] & BADPIX_STAR) {
	k = imageval[i] ;
	if (!flag_nointerpolate)output->image[i] = median[k] + sigimage[k]*(gasdev(&seed)) ; 
	/*	if (k ==0)printf ("%f %f %f\n",output->image[i],sigimage[k],median[k]); */ 
	/*	printf ("%d %d %f %f\n",k,i,median[k],sigimage[k]);	*/
      }   		   
    }
  /*temporarily turned bleeding trails off */
 
  if ((match > 0) && (exposure >10)) {
    if (horiztrails){
      sprintf(event,"Starting horizontal bleeding trails masking\n");
      reportevt(flag_verbose,STATUS,1,event);
      bleedtrailsx(output,centerx,centery,circlemean,radvalues,match,median,sigimage,1,flag_nointerpolate);
      bleedtrailsx(output,centerx,centery,circlemean,radvalues,match,median,sigimage,-1,flag_nointerpolate);}
    else
      {sprintf(event,"Starting vertical bleeding trails masking\n");
	reportevt(flag_verbose,STATUS,1,event);
	bleedingtrails(output,centerx,centery,circlemean,radvalues,match,1,flag_nointerpolate);
	bleedingtrails(output,centerx,centery,circlemean,radvalues,match,-1,flag_nointerpolate);}
  }
  free(tra);
  free(tdec);
  free(b2);
  free(r2);
  free(i2);
  free(circlemean);
  free(centerx);
  free(centery);
  free(radvalues);
  free(median);
  free(imageval);
}

float getmagnitude(band,B2MAG,R2MAG,NMAG)
     float B2MAG,R2MAG,NMAG;
     char  *band;
{
  float gmag,rmag,imag,zmag,Ymag,mag,magmax ;
  magmax=B2MAG;
  if (R2MAG>magmax)magmax=R2MAG;
  if (NMAG>magmax)magmax=NMAG;

  gmag = B2MAG - 0.15 ;
  rmag  =    R2MAG + 0.2 ;
  imag  = NMAG ;
  zmag =  imag - 0.025 - (0.8/1.9)*(rmag-imag) ;
  Ymag = zmag  + 0.22  - 0.408*(imag-zmag+0.5) ;
  switch ( *band ) {
  case 'g': mag = gmag ; break ;
  case 'r': mag = rmag ; break ;
  case 'i': mag = imag ; break ;
  case 'z': mag = zmag; break ;
  case 'Y': mag = Ymag; break ; 
  case 'U': mag = gmag; break ;
  case 'u': mag = gmag; break ;
  case 'un': mag = gmag; break ;
  case 'B': mag = gmag; break;
  case 'R':  mag = rmag; break;
  case 'V': mag = rmag; break ;
  case 'I': mag = imag; break;
  case 'It': mag = imag; break;
  case '337':  mag = gmag; break;
  case '390': mag = gmag; break ;
  case 'Bw':  mag = gmag; break ;
  case 'C' : mag = gmag; break ;
  case 'Gn' :  mag = gmag; break ;
  case '420' : mag = gmag; break ;
  case '454' : mag = gmag; break ;
  case 'o2' : mag = gmag; break ;
  case 'o3' : mag = gmag; break ;
  case 'O3' : mag = gmag; break ;
  case 'wh' : mag = gmag; break ;
  case 'wr475' : mag = gmag; break ;
  case 'wrc3' : mag = gmag; break ;
  case 'wrc4' : mag = gmag; break ;
  case 'wrhe2' : mag = gmag; break ;
  case 'Us' : mag = gmag; break ;
  case 'Ud' : mag = gmag; break ;
  case '493' : mag = rmag; break ;
  case '527' : mag = rmag; break ;
  case '579' : mag = rmag; break ;
  case '607' : mag = rmag; break ;
  case '666' :  mag = rmag; break ;
  case '705' : mag = rmag; break ;
  case 'D51' : mag = rmag; break ;
  case 'ha4' : mag = rmag; break ;
  case 'ha8' : mag = rmag; break ;
  case 'ha' : mag = rmag; break ;
  case 'ha12' : mag = rmag; break ;
  case 'ha16' : mag = rmag; break ;
  case 'M' : mag = rmag; break ;
  case 's2' : mag = rmag; break ;
  case 'VR' : mag = rmag; break ;
  case '755' : mag = imag; break ;
  case '802' : mag = imag; break ;
  case '815' : mag = imag; break ;
  case '823' : mag = imag; break ;
  case '848' : mag = imag; break ;
  case '918' : mag = imag; break ;
  case '973' : mag = imag; break ;

  default: mag=magmax;
  }
  return mag ;
}

void bleedingtrails(output,centerx,centery,circlemean,radvalues,match,step,flag_nointerpolate)
     desimage *output;
     int match,step,flag_nointerpolate ;
     float *radvalues,*centerx,*centery;
     double *circlemean;
{
  int xval,yval,trailmatch,i,j,k,trailx,traily,*masked,region[4],regionum,num,xmax,xmin,sum,count;
  float regionfactor,*regionsort,gasdev(),mode,sigma ;
  long seed = -30 ;
  void retrievescale(),pixelhisto() ;
  /* take care of bleeding trails */

  for (j=0;j<match;j++){
    num = 0 ;
    masked = (int *)calloc(20000,sizeof(int));
    traily = (int) (centery[j] + step * (radvalues[j] + 1)) ;
    trailx = (int) (centerx[j]) ;
    yval = traily ;
    trailmatch = YES ;
    xmax = centerx[j]; xmin = centerx[j] ;
    /*temporary fix for no stars in images corresponding to USNOB stars */
    if (circlemean[j] > 0)
      while ((yval <= output->axes[1])&& (yval>=0) && (trailmatch==YES) && 
	     (num<19000)){
	trailmatch = NO ; 
	for (xval=trailx-10;xval<trailx+10;xval++) {
	  i = (int) (yval*output->axes[0] + xval) ;
	  if (i>=0 && i<output->npixels)
	    {
	      if ((output->image[i] > circlemean[j]) ||
		  (fabs(output->image[i] - circlemean[j]) < sqrt(circlemean[j]))){
		if (xval>xmax)xmax = xval;
		if (xval<xmin)xmin =  xval ;
		output->mask[i] |= BADPIX_TRAIL ;
		output->varim[i] = 0 ;
		masked[num]=i ;

		num++ ;
		if (!flag_nointerpolate)output->image[i] = 0 ;
		trailmatch = YES ;
	      }
	    }
	}
	yval = yval + step ;
	/*      printf("%d\n",yval); */
      }
		
    if (num>0){
      region[0] = centerx[j]-3*(centerx[j]-xmin) ;
      region[1] = centerx[j]+3*(xmax - centerx[j]) ;
      region[2] = (int) (fminf((float) traily,(float) yval));
      region[3] = (int) (fmaxf((float) traily, (float) yval));

      if (region[0]<0)region[0]= 0;
      if (region[1]>output->axes[0])region[1]= output->axes[0];
      if (region[2]<0)region[2]= 0;
      if (region[3]>output->axes[1])region[3]= output->axes[1];


      regionum=(region[1]-region[0]+10)*
	(region[3]-region[2]+10);
      if (regionum>0){
	regionsort=(float *)calloc(regionum,sizeof(float));
	pixelhisto(output,region,regionsort,1,&regionfactor,&mode,&sigma);
	free(regionsort);
	if(fabs(regionfactor)<0.1)regionfactor = mode ;
	if (fabs(sigma)<0.1)sigma = sqrt(regionfactor);
	for (k=0;k<num;k++)
	  if(!flag_nointerpolate)output->image[masked[k]] = regionfactor + sigma*gasdev(&seed) ;
      }
			
    }
    free(masked);

  }
  return ;
}


void bleedtrailsx(output,centerx,centery,circlemean,radvalues,match,median,sigimage,step,flag_nointerpolate)
     desimage *output;
     int match,step,flag_nointerpolate;
     float *radvalues,*centerx,*centery,*sigimage,*median;
     double *circlemean;
{
  int xval,yval,trailmatch,i,j,k,trailx,traily,*masked,region[4],regionum,num,ymax,ymin;
  float regionfactor,*regionsort,gasdev(),mode,sigma ;
  long seed = -30 ;
  void retrievescale() ;
  /* take care of bleeding trails */

  for (j=0;j<match;j++){
    num = 0 ;
    masked = (int *)calloc(10000,sizeof(int));
    trailx = (int) (centerx[j] + step * (radvalues[j] + 1)) ;
    traily = (int) (centery[j]) ;
    xval = trailx ;
    trailmatch = YES ;
    ymax = centery[j]; ymin = centery[j] ;
    /*temporary fix for no stars in images corresponding to USNOB stars */
    if (circlemean[j] > 0)
      while ((xval < output->axes[0])&& (xval>=0) && (trailmatch==YES) && 
	     (num<9000)){
	trailmatch = NO ;               
	for (yval=traily-10;yval<traily+10;yval++) {
	  i = (int) (yval*output->axes[0] + xval) ;
	  if((i>=0) && (i<output->npixels)){
	    if ((output->image[i] > circlemean[j]) ||
		(fabs(output->image[i] - circlemean[j]) < sqrt(circlemean[j]))){
	      if (yval>ymax)ymax = yval;
	      if (yval<ymin)ymin =  yval ;
	      output->mask[i] |= BADPIX_TRAIL ;
	      output->varim[i] = 0 ;
	      masked[num]=i ;

	      num++ ;
	      output->image[i] = 0 ;
	      trailmatch = YES ;
	    }
	  }
	}
	xval = xval + step ;
	/*      printf("%d\n",yval); */
      }
    if (num>0){
      /*                        region[0] = centerx[j]-3*(centerx[j]-xmin) ;
				region[1] = centerx[j]+3*(xmax - centerx[j]) ;
				region[2] = (int) (fminf((float) traily,(float) yval));
				region[3] = (int) (fmaxf((float) traily, (float) yval));
      */
      region[2] = centery[j]-3*(centery[j]-ymin) ;                             
      region[3] = centery[j]+3*(ymax - centery[j]) ;                                                       
      region[0] = (int) (fminf((float) trailx,(float) xval));                                              
      region[1] = (int) (fmaxf((float) trailx, (float) xval));                                             
		     
  

      regionum=(region[1]-region[0])*
	(region[3]-region[2]);
      /*
	regionsort=(float *)calloc(regionum,sizeof(float));


	if (region[0]<0)region[0]= 0;
	if (region[0]>output->axes[0])region[0]= output->axes[0];
	if (region[1]<0)region[1]= 0;
	if (region[1]>output->axes[0])region[1]= output->axes[0];

	if (region[2]<0)region[2]= 0;
	if (region[3]>output->axes[1])region[3]= output->axes[1];
	retrievescale(output,region,regionsort,1,&regionfactor,&mode,&sigma);
	free(regionsort);
	if(fabs(regionfactor)<0.1)regionfactor = mode ;
	if (fabs(sigma)<0.1)sigma = sqrt(fabs(regionfactor));
	if (fabs(sigma)<0.1)sigma = sqrt(fabs(mode));
	printf("bleeding trails mask: median and mean =  %f %f",regionfactor,sigma); */
      free(regionsort);
      regionfactor = median[j] ;
      sigma = sigimage[j];
      for (k=0;k<num;k++)
	if(!flag_nointerpolate)output->image[masked[k]] = regionfactor + sigma*gasdev(&seed) ;
    }
    free(masked);
  }
  return ;
}




void  pixelhisto(image,scaleregionn,scalesort,flag_verbose,scalefactor,
		 mode,sigma)
     desimage *image;
     float	*scalesort,*scalefactor,*mode,*sigma;
     int	*scaleregionn,flag_verbose;
{
  int	i,x,y,loc,xdim,npixels,j,sum,count,regionum,
    jmax,jplus,jminus,num,histy[100];
  float	histx[100],ymax,fraction;
  void	shell(),reportevt();
  char	event[10000];

  i=0;
  xdim=image->axes[0];
  npixels=image->npixels;
  for (j=0;j<100;j++)
    {
      histy[j] = 0;
      histx[j] = 2*j - 99 ;
    } 
  count = 0;
  regionum=(scaleregionn[1]-scaleregionn[0])*
    (scaleregionn[3]-scaleregionn[2]);
  /* copy good image values into sorting vector */
  for (y=scaleregionn[2]-1;y<scaleregionn[3];y++) 
    for (x=scaleregionn[0]-1;x<scaleregionn[1];x++) {
      loc=x+y*xdim;
      if (loc>=0 && loc<npixels) {
	if (!(image->mask[loc])) { /* of pixel not masked */
	  if(i<regionum)scalesort[i]=image->image[loc];
	  i++;
	}
      }
    }
  if (i>=regionum)i=regionum-1 ;
  if (i<100) {
    sprintf(event,"Very small scale region = [%d:%d,%d:%d] & image= %s",
	    scaleregionn[0],scaleregionn[1],scaleregionn[2],scaleregionn[3],
	    image->name);
    reportevt(flag_verbose,QA,3,event);
    *scalefactor=0.0; /* mark image as unuseable */
  }
  else {
    /* sort list */
    shell(i,scalesort-1);
    /* grab the median */
    if (i%2) *scalefactor=scalesort[i/2];
    else *scalefactor=0.5*(scalesort[i/2]+scalesort[i/2-1]);
    if (scalesort[0]>100){
      for (num=0;num<i;num++)scalesort[num] -=  *scalefactor ;
    }
    /* build a histogram */
    for (num=0;num<i;num++)
      if(fabs(scalesort[num])<100)
	{
	  j= (scalesort[num]+100)/2;
	  histy[j] += 1 ;
	}

    for (j=0;j<100;j++)count+=histy[j];

    ymax=0.0;jmax=jplus=jminus=0;
    for (j=100;j>0;j--) {
      if (histy[j-1]>ymax) {
	ymax=histy[j-1];
	jmax=j-1;
      }	    
    }
    for (j=jmax;j<100;j++) 
      if (histy[j]<0.5*ymax) {
	jplus=j;
	break;
      }
    for (j=jmax;j>=0;j--) 
      if (histy[j]<0.5*ymax) {
	jminus=j;
	break;
      }

    *mode=histx[jmax];
    *sigma=histx[jplus]-histx[jminus];
    *sigma/=2.354;  /* change to sigma assuming gaussian distribution */

	      
    if (flag_verbose) {
      fraction= (float)i/(float)((scaleregionn[1]-scaleregionn[0]+1)*
				 (scaleregionn[3]-scaleregionn[2]+1));
      if (strncmp(image->name,"!",1)) sprintf(event,"Scale = [%d:%d,%d:%d] & fractional no of pixels= %.4f  & Scale=%.1f & Mode=%.1f & Sigma=%.2f & image=%s",
					      scaleregionn[0],scaleregionn[1],scaleregionn[2],
					      scaleregionn[3],fraction,*scalefactor,*mode,*sigma,image->name);
      else sprintf(event,"Scale =  [%d:%d,%d:%d] & fractional no of pixels = %.5f  & Scale=%.1f & Mode=%.1f & Sigma=%.2f & image=  %s",
		   scaleregionn[0],scaleregionn[1],scaleregionn[2],
		   scaleregionn[3],fraction,*scalefactor,*mode,*sigma,image->name+1);
      reportevt(flag_verbose,QA,1,event);
    }


  }
}

void  growcircleradius(image,radius,centerx,centery,skymean,circlemean,imageval,starno)
     desimage *image;
     float	*radius,centerx,centery,skymean,circlemean;
     int starno,*imageval;
     /* find distance between radius and corners of the image */
{
  float maxradius,th,krad,xmask,ymask,mean ;
  int i, ct, *badpixindex,num=0,j,xpos,ypos,*duplicatecheck; 
  maxradius = (*radius)*1.7;
  badpixindex=(int *)calloc(image->npixels,sizeof(int));
  duplicatecheck = (int *)calloc(image->npixels,sizeof(int));
  for (i=0;i<image->npixels;i++)duplicatecheck[i] = 0;
  ct = 0;
  for (krad = *radius;krad<maxradius;krad+=0.1)
    {
      mean = 0 ;
      num = 0 ;
 
      for(th= 0;th<=360;th+=0.25)
	{
	  xmask = krad*cos(th*3.1416/180.0);
	  ymask = krad*sin(th*3.1416/180.0); 	  		
	  xpos = (int) (xmask + centerx) ;
	  ypos = (int) (ymask + centery) ;
	  if(xpos<0)xpos=0.0 ;
	  if(ypos<0)ypos=0.0 ;
	  if(xpos >= image->axes[0])xpos=image->axes[0] -1 ;
	  if(ypos >= image->axes[1])ypos=image->axes[1] - 1;
	  i = ypos*image->axes[0] + xpos ;
	  if (!duplicatecheck[i])
	    {
	      duplicatecheck[i]= 1 ;	
	      badpixindex[ct]=i ;
	      ct++ ;
	      /* if (!image->mask[i]){num+= 1; mean += image->image[i]; printf("%f %d %f\n",mean,num,krad);} */
	      if (fabs(sin(th))<0.996)
		{num += 1; mean += image->image[i]; }
	      /* printf ("%d %d %d  %f %f %d %f\n",xpos,ypos,i,image->image[i],mean,num,sin(th)); */
	    }
	}
      mean /= num ;        
      /* printf ("%f %f %f %f\n",krad,mean,skymean,circlemean); */ 
	   
	   
      if ((mean<skymean) || (fabs(mean-skymean)<sqrt(skymean)))break ;
    }


  if ((krad-(*radius)) > 0.4)
    {
      for (j=0;j<ct;j++)
	{ 
	  if (!(image->mask[badpixindex[j]] & BADPIX_STAR))
	    {
	      image->mask[badpixindex[j]] |= BADPIX_STAR;
	      imageval[badpixindex[j]] = starno ;
	    }
	}
    }
  *radius = krad ;

  free(badpixindex);
  free(duplicatecheck);
  return ;
}								       
