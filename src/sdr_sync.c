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
#include "input_sdr.h" // for samples_fft_swap
#include "tunables.h"

#define dbg 0

double mag_squared(fftw_complex sample) {
    double x = sample[0];
    double y = sample[1];
    return x * x + y * y;
}

double cpx_abs(fftw_complex sample) {
  return sqrt(mag_squared(sample));
}

void cpx_mul(fftw_complex a, fftw_complex b, float bConjSign, fftw_complex *res) {
  (*res)[0] = a[0]*b[0]           - a[1]*bConjSign*b[1];
  (*res)[1] = a[0]*bConjSign*b[1] + a[1]*b[0];
}


// another aproach to time sync:
// step > 1 tries to make things faster, but somewhat less reliable.
// if sync is good already, it's somewhat less robust than dab_fine_time_sync,
// even with step=1.


// start can be negative up to >=-DAB_T_NULL
int time_sync_sub(fftw_complex * frame, double * rmin, double * rmax, int start, int end, int step)
{
  // assume most DC offset is removed by hardware
  double sum = 0.0;
  int i;

  if ( rmin != NULL || rmax != NULL) {
    for (i = start; i < end; i += step) {
      int il = (i + DAB_T_FRAME) % DAB_T_FRAME;
      sum += mag_squared(frame[il]);
    }
  }

  int minPos = start;
  double minVal = sum;
  double maxVal = sum;
  for (i = start; i < end; i += step) {
    int il = (i + DAB_T_FRAME) % DAB_T_FRAME;
    int ir = (i + DAB_T_NULL) % DAB_T_FRAME;
    sum += mag_squared(frame[ir]) - mag_squared(frame[il]);
    if (sum < minVal) {
      /* fprintf(stderr, ";%6d-%d ", minPos, i-minPos); // old minPos */
      minVal = sum;
      minPos = i;
    }
    else if (sum > maxVal) {
      maxVal = sum;
    }
  }

  if (rmin != NULL) { *rmin = minVal; }
  if (rmax != NULL) { *rmax = maxVal; }
  /* fprintf(stderr, "SUB MINPOS/VAL : %6d %8.3f %8.3f\n", 2*minPos, minVal, maxVal); */
  return minPos;
}

int time_sync_still_good_2(fftw_complex * frame) {
  double min, max;
  int ts = time_sync_sub(frame, &min, &max, -DAB_T_NULL/2, DAB_T_NULL/2, 1);
  return abs(ts) < DAB_T_NULL/4 && min < REALATIVE_ENERGY_THRESHOLD * max;
}

int time_sync_full(fftw_complex * frame) 
{
  int step = TIME_SYNC_INIT_STEP;
  int neigh = 3;
  int ts1 = time_sync_sub(frame, NULL, NULL, 0, DAB_T_FRAME - neigh*step, step);
  if (ts1 > DAB_T_FRAME - DAB_T_NULL + neigh*step) {
    ts1 -= DAB_T_FRAME;
  }
  int ts2 = time_sync_sub(frame, NULL, NULL, ts1 - neigh*step, ts1 + neigh*step, 1);
  return 2*ts2; // bytes
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


void debug_time_sync(fftw_complex * frame)
{
  int i;
  fprintf(stderr, "decr steps forward: ");
  float suml = 0;
  float sumr = 0;
  int lastprinted = 0;
  for (i=0; i<DAB_T_SYM; i++) {
    suml += cpx_abs(frame[i]);
    sumr += cpx_abs(frame[DAB_T_NULL + i]);
    if (sumr < suml) {
      if (!lastprinted) {
	fprintf(stderr, "%d", i);
	lastprinted = 1;
      }
    }
    else {
      if (lastprinted) {
	fprintf(stderr, "..%d, ", i);
	lastprinted = 0;
      }
    }
  }
  fprintf(stderr, "\n");

  fprintf(stderr, "decr steps backward: ");
  suml = 0;
  sumr = 0;
  lastprinted = 0;
  for (i=1; i<DAB_T_SYM; i++) {
    suml += cpx_abs(frame[DAB_T_FRAME - i]);
    sumr += cpx_abs(frame[DAB_T_NULL - i]);
    if (suml < sumr) {
      if (!lastprinted) {
	fprintf(stderr, "-%d", i);
	lastprinted = 1;
      }
    }
    else {
      if (lastprinted) {
	fprintf(stderr, "..-%d, ", i);
	lastprinted = 0;
      }
    }
  }
  fprintf(stderr, "\n");
}



// determines fine time shift (within symbol) by convoluting the PRS
// frame: start of frame with (previous) timeshifts applied -> starts with null symbol
// ret: new timeshift (in bytes to skip = 2*samples) relative to the frame input
int dab_fine_time_sync(fftw_complex * frame)
{
  // debug_time_sync(frame);

  // We want to shift backwards a bit so we're always inside the symbol.
  int guardShift = DAB_T_GUARD/2;
  // guard (with respedt to PRS) is considered a the beginning of the
  // symbol (TODO: no reference found).  
  int start = DAB_T_NULL + DAB_T_GUARD - guardShift;
  /* correlation in frequency domain 
     e.g. J.Cho "PC-based receiver for Eureka-147" 2001
     e.g. K.Taura "A DAB receiver" 1996
  */
  debug_print_array(DAB_T_SYM, frame + start, "prs_received.dat");

  /* first we have to transfer the receive prs symbol in frequency domain */
  fftw_complex prs_received_fft[DAB_T_CS];
  samples_fft_swap(&frame[start], prs_received_fft);
  debug_print_array(DAB_T_CS, prs_received_fft, "prs_received_fft.dat");

  /* fftshift the received prs
     at this point we have to be coarse frequency sync 
     however we can simply shift the bins */
  fftw_complex prs_rec_shift[DAB_CARRIERS];
  int i;
  int kk=0;
  for (i=0; i<DAB_T_CS; i++) {
    int empty_channels = DAB_T_CS - DAB_CARRIERS;
    // 256..1792 \ 1024
    if (empty_channels/2 <= i && i != DAB_T_CS/2 && i <= DAB_CARRIERS+empty_channels/2) {
      // no frequency deinterlacing here
      prs_rec_shift[kk][0] = prs_received_fft[i][0];
      prs_rec_shift[kk][1] = prs_received_fft[i][1];
      kk++;
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

  int maxPos=0;
  float maxVal=-99999;
  for (i=0; i<DAB_CARRIERS; i++) {
    float tempVal = cpx_abs(convoluted_prs_time[i]);
    if (tempVal>maxVal) {
      maxPos = i;
      maxVal = tempVal;
    }
  }

  int shiftSamples = maxPos;
  // rescale to full DAB_T_CS
  shiftSamples = (shiftSamples * DAB_T_CS) / DAB_CARRIERS;
  // undo the guardShift
  shiftSamples = (shiftSamples - guardShift) % DAB_T_CS;
  // wrap around to +/- DAB_T_CS
  if (shiftSamples >= DAB_T_CS/2) {
    shiftSamples -= DAB_T_CS;
  }

  // TODO: for unknown reason this gives a shift of -34..-30 bytes every time

  return shiftSamples * 2; // convert to bytes to skip
}




// idea: find center of mas of active frequencies, avarage over all frames
// return frequency shift in DAB_F_C(=kHz)
int dab_coarse_freq_sync_4(int symcount, fftw_complex (*symbols)[DAB_T_CS])
{
  int usedWidth = DAB_CARRIERS + 1;  // the middle frequency (tuner mixing) is not used
  int middlePos = DAB_CARRIERS/2;  // i.e.: somewhat towards lower indices

  double sum = 0;
  double maxPos = middlePos;
  double maxVal = sum;

  int i;
  for (i = 0; i < DAB_T_CS - usedWidth; i++) {
    int isym;
    for (isym = 0; isym < symcount; isym++) {
      double vl = cpx_abs(symbols[isym][i]);
      double vr = cpx_abs(symbols[isym][i   + usedWidth]);
      double dr1 = vr - vl;
    
      sum += dr1;
    }

    // middle pos is usually noise (from mixing), don't catch that:
    // removed - cpx_abs(syms[i + middlePos]);
    double val = sum;
    double pos = i + middlePos + 1 - DAB_T_CS/2; // +1: steps above prepare i+1 from i
    
    if (val > maxVal) {
      maxPos = pos;
      maxVal = val;
    }
  }

  return maxPos;
}


// determines phase shift over one or more symbols; see (guess) https://en.wikipedia.org/wiki/Cyclic_prefix
// ret: frequency shift in Hz, domain is +/-500
double dab_fine_freq_corr(fftw_complex * dab_frame)
{
  int numSymbols = 15;
  fftw_complex sum = {0.0, 0.0};
  double angle, ffs;
  int i, j;

  for (j=DAB_SYMBOLS_IN_FRAME-numSymbols; j<DAB_SYMBOLS_IN_FRAME; j++) {
    // a frequency shift might have been applied in the middle of the frame
    // determine the new shift from the newest data (last symbols)
    int skip = DAB_T_NULL +  j * DAB_T_SYM;
    for (i=0; i<DAB_T_GUARD; i++) {
      fftw_complex lr;    
      // the guard is cyclic, so we can calculate the phase shift between them
      cpx_mul(dab_frame[skip+DAB_T_CS+i], dab_frame[skip+i], -1, &lr);
      sum[0] += lr[0];
      sum[1] += lr[1];
    }
  }

  angle = atan2(sum[1], sum[0]);

  ffs = angle / (2 * M_PI) * DAB_F_C;
  /* fprintf(stderr, "\n%f\n",ffs); */

  return ffs;
}
