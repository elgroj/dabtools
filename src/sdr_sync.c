/*
This file is part of rtl-dab
trl-dab is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Foobar is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with rtl-dab.  If not, see <http://www.gnu.org/licenses/>.


david may 2012
david.may.muc@googlemail.com

*/

#include "sdr_sync.h"
#include "sdr_prstab.c"
#include "dab_constants.h"

#define dbg 0

float mag_squared(fftw_complex sample) {
    float x = sample[0];
    float y = sample[1];
    return x * x + y *y;
}


// see also (maybe?) https://en.wikipedia.org/wiki/Guard_interval
// real: real part
// filt: temporary buffer to do work
// force_timesync: boolean, timesync is done anyway; (why? where to?)
// ret: amount (bytes?) of shift; should be multiple of 20
uint32_t dab_coarse_time_sync(int8_t * real, float * filt, uint8_t force_timesync) {
  int32_t tnull = DAB_T_NULL; // was 2662? why?   (2656)
  int32_t j,k;

  // check for energy in fist tnull samples
  float e=0;
  float threshold=5000;
  for (k=0; k<tnull; k+=10) {
    e += (float) abs(real[k]);
  }
  //fprintf(stderr, "Energy over nullsymbol: %f\n", e);
#if dbg
  fprintf(stderr,"Energy over nullsymbol: %f\n",e);
#endif
  if (e<threshold && (force_timesync==0))
    return 0;
  //fprintf(stderr,"Resync\n");
  // energy was to high so we assume we are not in sync
  // subsampled filter to detect where the null symbol is
  for (j=0; j<(DAB_T_FRAME-tnull)/10; j++)
    filt[j] = 0;
  for (j=0; j<DAB_T_FRAME-tnull; j+=10)
    for (k=0;k<tnull;k+=10)
      filt[j/10] += (float) abs(real[j+k]);

  // finding the minimum in filtered data gives position of null symbol
  float minVal=9999999;
  float maxVal=0.0;//
  uint32_t minPos=0;
  for (j=0; j<(DAB_T_FRAME-tnull)/10; j++) {
    if (filt[j]<minVal) {
      minVal = filt[j];
      minPos = j*10;
    }
    // debug
    if (filt[j]>maxVal) {
      maxVal = filt[j];
    }
  }
  fprintf(stderr, "minVal %7.2f, maxVal %7.2f, minPos*2 %4d\n", minVal, maxVal, minPos*2);

  //fprintf(stderr,"calculated position of nullsymbol: %f",minPos*2);
  return minPos*2;
}


int32_t dab_fine_time_sync(fftw_complex * frame){
  /* correlation in frequency domain 
     e.g. J.Cho "PC-based receiver for Eureka-147" 2001
     e.g. K.Taura "A DAB receiver" 1996
  */
#if dbg
  int k;
  FILE *fh0;
  fh0 = fopen("prs_received.dat","w+");
  for(k=0; k<DAB_T_SYM; k++) {
    fprintf(fh0,"%f\n",(frame[DAB_T_NULL+k][0]));
    fprintf(fh0,"%f\n",(frame[DAB_T_NULL+k][1]));
  }
  fclose(fh0);
#endif

  /* first we have to transfer the receive prs symbol in frequency domain */
  fftw_complex prs_received_fft[DAB_T_CS];
  fftw_plan p;
  p = fftw_plan_dft_1d(DAB_T_CS, &frame[DAB_T_NULL+DAB_T_GUARD], &prs_received_fft[0], 
		       FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
#if dbg
  FILE *fh1;
  fh1 = fopen("prs_received_fft.dat","w+");
  for(k=0; k<DAB_T_CS; k++) {
    fprintf(fh1,"%f\n",(prs_received_fft[k][0]));
    fprintf(fh1,"%f\n",(prs_received_fft[k][1]));
  }
  fclose(fh1);
#endif
  /* now we build the complex conjugate of the known prs */
  // 1536 as only the carries are used = DAB_CARRIERS
  fftw_complex prs_star[DAB_CARRIERS];
  int i;
  for (i=0; i<DAB_CARRIERS; i++) {
    prs_star[i][0] = prs_static[i][0];
    prs_star[i][1] = -1 *  prs_static[i][1];
  }
#if dbg
  FILE *fh2;
  fh2 = fopen("prs_star.dat","w+");
  for(k=0; k<DAB_CARRIERS; k++) {
    fprintf(fh2,"%f\n",(prs_star[k][0]));
    fprintf(fh2,"%f\n",(prs_star[k][1]));
  }
  fclose(fh2);
#endif


  /* fftshift the received prs
     at this point we have to be coarse frequency sync 
     however we can simply shift the bins */
  fftw_complex prs_rec_shift[DAB_CARRIERS];
  // TODO allow for coarse frequency shift !=0 
  // int32_t cf_shift = 0;
  // matlab notation (!!!-1)
  // 769:1536+s   DAB_CARRIERS+s
  //  2:769+s why 2? I dont remember, but peak is very strong
  int mid = DAB_CARRIERS/2; // 768
  for (i=0; i<DAB_CARRIERS; i++) {
    if (i < mid) {
      prs_rec_shift[i][0] = prs_received_fft[i+1280][0];
      prs_rec_shift[i][1] = prs_received_fft[i+1280][1];
    }
    else { // i>=768
      prs_rec_shift[i][0] = prs_received_fft[i-765][0];
      prs_rec_shift[i][1] = prs_received_fft[i-765][1];
    }
  }
#if dbg
  FILE *fh3;
  fh3 = fopen("prs_rec_shift.dat","w+");
  for(k=0; k<DAB_CARRIERS; k++) {
    fprintf(fh3,"%f\n",(prs_rec_shift[k][0]));
    fprintf(fh3,"%f\n",(prs_rec_shift[k][1]));
  }
  fclose(fh3);
#endif


  
  /* now we convolute both symbols */
  fftw_complex convoluted_prs[DAB_CARRIERS];
  int s;
  for (s=0; s<DAB_CARRIERS; s++) {
    convoluted_prs[s][0] = prs_rec_shift[s][0]*prs_star[s][0] - prs_rec_shift[s][1]*prs_star[s][1];
    convoluted_prs[s][1] = prs_rec_shift[s][0]*prs_star[s][1] + prs_rec_shift[s][1]*prs_star[s][0];
  }

  /* and finally we transfer the convolution back into time domain */
  fftw_complex convoluted_prs_time[DAB_CARRIERS]; 
  fftw_plan px;
  px = fftw_plan_dft_1d(DAB_CARRIERS, &convoluted_prs[0], &convoluted_prs_time[0], 
			FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(px);
  fftw_destroy_plan(px);
#if dbg
  FILE *fh4;
  fh4 = fopen("convoluted_prs_time.dat","w+");
  for(k=0; k<DAB_CARRIERS; k++) {
    fprintf(fh4,"%f\n",(convoluted_prs_time[k][0]));
    fprintf(fh4,"%f\n",(convoluted_prs_time[k][1]));
  }
  fclose(fh4);
#endif

  uint32_t maxPos=0;
  float tempVal = 0;
  float maxVal=-99999;
  for (i=0; i<DAB_CARRIERS; i++) {
    tempVal = sqrt((convoluted_prs_time[i][0]*convoluted_prs_time[i][0])
		   + (convoluted_prs_time[i][1]*convoluted_prs_time[i][1]));
    if (tempVal>maxVal) {
      maxPos = i;
      maxVal = tempVal;
    }
  }

#if dbg
  fprintf(stderr,"Fine time shift: %d\n",maxPos);
#endif
  if (maxPos<DAB_CARRIERS/2) {
    return maxPos*2+16;
  } else {
    return ((maxPos-(DAB_CARRIERS))*2);
  }
  //return 0;
}





// return -14..14 freq shift in kHz
int32_t dab_coarse_freq_sync_2(fftw_complex * symbols){
  int len = 128;
  fftw_complex convoluted_prs[len];
  int s;
  int freq_hub = 14; // + and - center freq
  int k;
  float global_max = -99999;
  int global_max_pos=0; 
  for (k=-freq_hub;k<=freq_hub;k++) {
    
    // complex element-wise product over len symbol-frequencies with offset 256
    // against some static frequency pattern (sdr_prstab.c, the subband-structure?)
    // (->gives convolution when fft-ed again)
    // freq_hub+k goes 0..2*freq_hub
    for (s=0;s<len;s++) {
      convoluted_prs[s][0] = prs_static[freq_hub+s][0]*symbols[freq_hub+k+256+s][0]-
	(-1)*prs_static[freq_hub+s][1]*symbols[freq_hub+k+256+s][1];
      convoluted_prs[s][1] = prs_static[freq_hub+s][0]*symbols[freq_hub+k+256+s][1]+
	(-1)*prs_static[freq_hub+s][1]*symbols[freq_hub+k+256+s][0];
    }
    fftw_complex convoluted_prs_time[len]; 
    fftw_plan px;
    px = fftw_plan_dft_1d(len, &convoluted_prs[0], &convoluted_prs_time[0], FFTW_BACKWARD, FFTW_ESTIMATE);
    fftw_execute(px);
    fftw_destroy_plan(px);
    
    
#if dbg
    FILE *fh0;
    fh0 = fopen("convoluted_prs_coarse.dat","w+");
    for(s=0;s<len;s++) {
      fprintf(fh0,"%f\n",(convoluted_prs_time[s][0]));
      fprintf(fh0,"%f\n",(convoluted_prs_time[s][1]));
    }
    fclose(fh0);
#endif
    
    // uint32_t maxPos=0;
    float tempVal = 0;
    float maxVal=-99999;
    for (s=0;s<len;s++) {
      // complex abs value
      tempVal = sqrt((convoluted_prs_time[s][0]*convoluted_prs_time[s][0])+(convoluted_prs_time[s][1]*convoluted_prs_time[s][1]));
      // maxpos is the *time*-shift where correlation is best
      // also, it is not used :-)
      if (tempVal>maxVal) {
	// maxPos = s;
	maxVal = tempVal;
      }
    }
    //fprintf(stderr,"%f ",maxVal);
    
    
    if (maxVal>global_max) {
      global_max = maxVal;
      global_max_pos = k;
    }
  }
  //fprintf(stderr,"MAXPOS %d\n",global_max_pos);
  return global_max_pos;
}


// see (guess) https://en.wikipedia.org/wiki/Cyclic_prefix
double dab_fine_freq_corr(fftw_complex * dab_frame,int32_t fine_timeshift){
  fftw_complex *left; // = early in time
  fftw_complex *right; // = arrived somewhat later
  fftw_complex *lr;
  double angle[DAB_T_GUARD];
  double mean=0;
  double ffs;
  left = fftw_malloc(sizeof(fftw_complex) * DAB_T_GUARD);
  right = fftw_malloc(sizeof(fftw_complex) * DAB_T_GUARD);
  lr = fftw_malloc(sizeof(fftw_complex) * DAB_T_GUARD);
  uint32_t i;
  fine_timeshift = 0;  // TODO why are we overriding a value argument???
  for (i=0; i<DAB_T_GUARD; i++) {
    left[i][0]  = dab_frame[DAB_T_NULL+DAB_T_CS+i+fine_timeshift][0];
    left[i][1]  = dab_frame[DAB_T_NULL+DAB_T_CS+i+fine_timeshift][1];
    right[i][0] = dab_frame[DAB_T_NULL         +i+fine_timeshift][0];
    right[i][1] = dab_frame[DAB_T_NULL         +i+fine_timeshift][1];
  }
  // (complex) element-wise product; phase of lr is sum of phases
  for (i=0; i<DAB_T_GUARD; i++) {
    lr[i][0] = left[i][0]*right[i][0] - left[i][1]*(-1)*right[i][1];
    lr[i][1] = left[i][0]*(-1)*right[i][1] + left[i][1]*right[i][0];
  }
  
  for (i=0; i<DAB_T_GUARD; i++) {
   angle[i] = atan2(lr[i][1],lr[i][0]);
  }
  for (i=0; i<DAB_T_GUARD; i++) {
    mean = mean + angle[i];
  }
  mean = (mean/DAB_T_GUARD);
  //printf("\n%f %f\n",left[0][0],left[0][1]);
  //printf("\n%f %f\n",right[0][0],right[0][1]);
  //printf("\n%f %f\n",lr[0][0],lr[0][1]);
  //printf("\n%f\n",angle[0]);

  ffs = mean / (2 * M_PI) * 1000;
  //printf("\n%f\n",ffs);

  fftw_free(left);
  fftw_free(right);
  fftw_free(lr);
    
 return ffs;
}
