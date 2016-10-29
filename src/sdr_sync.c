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
#include "input_sdr.h" // for FILT_SAMPLING

#define dbg 0

float mag_squared(fftw_complex sample) {
    float x = sample[0];
    float y = sample[1];
    return x * x + y * y;
}

float cpx_abs(fftw_complex sample) {
  return sqrt(mag_squared(sample));
}

void cpx_mul(fftw_complex a, fftw_complex b, float bConjSign, fftw_complex *res) {
  (*res)[0] = a[0]*b[0]           - a[1]*bConjSign*b[1];
  (*res)[1] = a[0]*bConjSign*b[1] + a[1]*b[0];
}


// see also (maybe?) https://en.wikipedia.org/wiki/Guard_interval
// real: real part
// filt: temporary buffer to do work
// force_timesync: boolean, timesync is done anyway; (why? where to?)
// ret: bytes (2*samples) of shift; should be multiple of 20
uint32_t dab_coarse_time_sync(int8_t * real, float * filt, uint8_t force_timesync) {
  int32_t tnull = DAB_T_NULL; // was 2662? why?   (2656)
  int32_t j,k;

  // check for energy in fist tnull samples
  float e=0;
  for (k=0; k<tnull; k+=FILT_SAMPLING) {
    e += (float) abs(real[k]);
  }

  // collect same amount of hashed (kind of random) samples for comparison
  float eFullFrame = 0;
  j=3331;
  for (k=0; k<tnull; k+=FILT_SAMPLING) {
    eFullFrame += (float) abs(real[j % DAB_T_FRAME]);
    j = (j + 31)*331 % 33333331;  // note: 33*1 used here are all prime
  }

  float ratio = e/eFullFrame;
#if dbg
  fprintf(stderr,"Energy over nullsymbol: %f  full: %f  ratio: %f\n", e, eFullFrame, ratio);
#endif
  /* float threshold=5000; // FIXME: fixed threshold, this looks stupid to me */
  /* if (e<threshold && (force_timesync==0)) */
  /*   return 0; */
  float realativeThreshold = 0.3;
  if (ratio < realativeThreshold && !force_timesync) {
    return 0;
  }
  fprintf(stderr,"Energy over nullsymbol: %f  full: %f  ratio: %f\n", e, eFullFrame, ratio);


  //fprintf(stderr,"Resync\n");
  // energy was to high so we assume we are not in sync
  // subsampled filter to detect where the null symbol is
  for (j=0; j<(DAB_T_FRAME-tnull)/FILT_SAMPLING; j++) {
    filt[j] = 0;
  }
  int tnullSamp; // ((tnull-1)/FILT_SAMPLING)*FILT_SAMPLING;
  for (j=0; j<tnull; j+=FILT_SAMPLING) {
    filt[0] += (float) abs(real[j]);
    tnullSamp = j; // easy way: just record last j
  }
  for (j=FILT_SAMPLING; j<DAB_T_FRAME-tnull; j+=FILT_SAMPLING) {
    int fi = j/FILT_SAMPLING;  // starts with 1
    filt[fi] = filt[fi-1] 
      - (float) abs(real[j - FILT_SAMPLING])
      + (float) abs(real[j + tnullSamp]);
  }
  /* for (j=0; j<DAB_T_FRAME-tnull; j+=FILT_SAMPLING) */
  /*   for (k=0; k<tnull; k+=FILT_SAMPLING) */
  /*     filt[j/FILT_SAMPLING] += (float) abs(real[j+k]); */

  // finding the minimum in filtered data gives position of null symbol
  float minVal=9999999;
  float maxVal=0.0;//
  uint32_t minPos=0;
  for (j=0; j<(DAB_T_FRAME-tnull)/FILT_SAMPLING; j++) {
    if (filt[j]<minVal) {
      minVal = filt[j];
      minPos = j*FILT_SAMPLING;
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

void debug_print_array(int len, fftw_complex *arr, char* fileName)
{
  int i;
  FILE *fh0;
  if (dbg) {
    fh0 = fopen(fileName, "w+");
    for(i=0; i<len; i++) {
      fprintf(fh0,"%f\n",(arr[i][0]));
      fprintf(fh0,"%f\n",(arr[i][1]));
    }
    fclose(fh0);
  }
}


// frame: start of frame with (previous) timeshifts applied -> starts with null symbol
// ret: new timeshift (in bytes = 2*samples) relative to the frame input
int32_t dab_fine_time_sync(fftw_complex * frame) 
{
  /* correlation in frequency domain 
     e.g. J.Cho "PC-based receiver for Eureka-147" 2001
     e.g. K.Taura "A DAB receiver" 1996
  */
  debug_print_array(DAB_T_SYM-DAB_T_NULL, frame + DAB_T_NULL, "prs_received.dat");

  /* first we have to transfer the receive prs symbol in frequency domain */
  fftw_complex prs_received_fft[DAB_T_CS];
  fftw_plan p;
  p = fftw_plan_dft_1d(DAB_T_CS, &frame[DAB_T_NULL+DAB_T_GUARD], &prs_received_fft[0], 
		       FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  debug_print_array(DAB_T_CS, prs_received_fft, "prs_received_fft.dat");

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
  int i;
  for (i=0; i<DAB_CARRIERS; i++) {
    // 1280 = 2048-768 = DAB_T_CS - mid
    if (i < mid) {
      // prs_received_fft[1280..2048-1] -> prs_rec_shift[0..768-1]
      prs_rec_shift[i][0] = prs_received_fft[i+(DAB_T_CS-mid)][0];
      prs_rec_shift[i][1] = prs_received_fft[i+(DAB_T_CS-mid)][1];
    }
    else { // i>=768
      // FIXME: magic numbers here: 765 = 768-3
      int unexplained_shift = 3; // was 3
      // prs_received_fft[3..771-1] -> prs_rec_shift[768..1536-1]
      prs_rec_shift[i][0] = prs_received_fft[i-(mid-unexplained_shift)][0];
      prs_rec_shift[i][1] = prs_received_fft[i-(mid-unexplained_shift)][1];
    }
  }
  debug_print_array(DAB_CARRIERS, prs_rec_shift, "prs_rec_shift.dat");

  
  /* now we convolute both symbols */
  fftw_complex convoluted_prs[DAB_CARRIERS];
  int s;
  for (s=0; s<DAB_CARRIERS; s++) {
    cpx_mul(prs_rec_shift[s], prs_static[s], -1, &convoluted_prs[s]);
  }

  /* and finally we transfer the convolution back into time domain */
  fftw_complex convoluted_prs_time[DAB_CARRIERS]; 
  fftw_plan px;
  px = fftw_plan_dft_1d(DAB_CARRIERS, &convoluted_prs[0], &convoluted_prs_time[0], 
			FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute(px);
  fftw_destroy_plan(px);
  debug_print_array(DAB_CARRIERS, convoluted_prs_time, "convoluted_prs_time.dat");

  uint32_t maxPos=0;
  float tempVal = 0;
  float maxVal=-99999;
  for (i=0; i<DAB_CARRIERS; i++) {
    tempVal = cpx_abs(convoluted_prs_time[i]);
    if (tempVal>maxVal) {
      maxPos = i;
      maxVal = tempVal;
    }
  }

  // FIXME: magic number here: 8 (was 16 before refactoring of *2)
  int magicAdditiveOnlyPosSide = 8;
  int shiftSamples =
    (maxPos < DAB_CARRIERS/2)
    ? maxPos + magicAdditiveOnlyPosSide
    : maxPos - DAB_CARRIERS;

  /* int oldMethod =  dab_fine_time_sync_old(frame); */
  /* fprintf(stderr, "FINE TIME SYNC COMPARISON: %10d %10d %10d \n",  */
  /* 	  oldMethod, shiftSamples*2, oldMethod-shiftSamples*2); */

  return shiftSamples * 2;
}





// return -14..14 freq shift in kHz
int32_t dab_coarse_freq_sync_2(fftw_complex * symbols)
{
  // FIXME: magic numbers here: 128, 14, also 256
  int len = 128;
  fftw_complex convoluted_prs[len];
  int s;
  int freq_hub = 14; // + and - center freq
  int k;
  float global_max = -99999;
  int global_max_pos=0; 

  for (k=-freq_hub;k<=freq_hub;k++) {
    // (->gives convolution when fft-ed again)
    // freq_hub+k goes 0..2*freq_hub
    for (s=0;s<len;s++) {
      cpx_mul(symbols[freq_hub+k+256+s], prs_static[freq_hub+s], -1, &convoluted_prs[s]);
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
      tempVal = cpx_abs(convoluted_prs_time[s]);
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
// ret: frequency shift in Hz
double dab_fine_freq_corr(fftw_complex * dab_frame)
{
  fftw_complex *left;
  fftw_complex *right;
  fftw_complex *lr;
  // double angle[DAB_T_GUARD];
  double mean = 0;
  double meanofsq =0;
  double ffs;
  left = fftw_malloc(sizeof(fftw_complex) * DAB_T_GUARD);
  right = fftw_malloc(sizeof(fftw_complex) * DAB_T_GUARD);
  lr = fftw_malloc(sizeof(fftw_complex) * DAB_T_GUARD);
  uint32_t i;
  for (i=0; i<DAB_T_GUARD; i++) {
    // fine_timeshift+DAB_T_NULL: we are looking at the Phase Reference Symbol
    left[i][0]  = dab_frame[DAB_T_NULL+DAB_T_CS+i][0];
    left[i][1]  = dab_frame[DAB_T_NULL+DAB_T_CS+i][1];
    right[i][0] = dab_frame[DAB_T_NULL         +i][0];
    right[i][1] = dab_frame[DAB_T_NULL         +i][1];
  }

  // the guard is cyclic, so we can calculate the phase shift between them
  for (i=0; i<DAB_T_GUARD; i++) {
    cpx_mul(left[i], right[i], -1, &lr[i]);
  }
  
  for (i=0; i<DAB_T_GUARD; i++) {
    // angle[i] = atan2(lr[i][1], lr[i][0]);
    double angle = atan2(lr[i][1], lr[i][0]);
    mean += angle;
    meanofsq += angle*angle;
  }
  mean /= DAB_T_GUARD;
  meanofsq /= DAB_T_GUARD;

  // magic numbers here: 1000  (samp_rate / DAB_T_CS)
  ffs = mean / (2 * M_PI) * 1000;
  // fprintf(stderr, "mean phase: %6.3f  +/-%7.3f;  ffs: %4.1f\n", 
  //    mean, sqrt(meanofsq - (mean*mean)), ffs);
  //fprintf(stderr, "\n%f\n",ffs);

  fftw_free(left);
  fftw_free(right);
  fftw_free(lr);
    
  return ffs;
}
