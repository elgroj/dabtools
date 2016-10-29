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

#include <assert.h>

#include "dab.h"
#include "dab_tables.h"
#include "input_sdr.h"
#include "sdr_sync.h"


void swap_blocks(int len, fftw_complex* sym) 
{
  fftw_complex tmp;
  int i;
  for (i = 0; i < len/2; i++) {
    tmp[0] = sym[i][0];
    tmp[1] = sym[i][1];
    sym[i][0] = sym[i + len/2][0];
    sym[i][1] = sym[i + len/2][1];
    sym[i + len/2][0] = tmp[0];
    sym[i + len/2][1] = tmp[1];
  }
}


void dqpsk_step(fftw_complex a, fftw_complex b, fftw_complex *out) 
{
  double magsqb = b[0]*b[0] + b[1]*b[1];
  (*out)[0] = (a[0]*b[0] + a[1]*b[1]) / magsqb;
  (*out)[1] = (a[0]*b[1] - a[1]*b[0]) / magsqb;
}

int sdr_demod(struct demapped_transmission_frame_t *tf, struct sdr_state_t *sdr) {
  int i,j;

  assert(DAB_T_FRAME == DAB_SYMBOLS_IN_FRAME * DAB_T_SYM + DAB_T_NULL);
  assert(DAB_T_SYM == DAB_T_CS + DAB_T_GUARD);

  //fprintf(stderr, "sdr_demod\n");
  //fprintf(stderr, ".");
  tf->has_fic = 0;

  /* resetting coarse freqshift */
  sdr->coarse_freq_shift = 0;
  
  /* write input data into fifo */
  for (i=0;i<sdr->input_buffer_len;i++) {
    cbWrite(&(sdr->fifo),&sdr->input_buffer[i]);
  }

  // timeshifts in bytes
  int offset = sdr->coarse_timeshift + sdr->fine_timeshift;

  /* Check for data in fifo */
  if ((int)sdr->fifo.count < 2*DAB_T_FRAME + offset) {
    // fprintf(stderr, "   no data in fifo %d < %d\n", sdr->fifo.count, DAB_T_FRAME * 3);
    return 0;
  }

  fprintf(stderr, "offset: %6d; cts: %6d fts: %6d\n",
  	  offset, sdr->coarse_timeshift, sdr->fine_timeshift);
  
  /* read fifo, 2 reads (Q/I) per sample */
  sdr_read_fifo(&(sdr->fifo), 2*DAB_T_FRAME, offset, sdr->buffer);
  sdr->fine_timeshift = 0; // time shift was part of offset

  
  /* give the AGC some time to settle */
  if (sdr->startup_delay<=GAIN_SETTLE_TIME) {
    // (frame has been read and will be thrown away)
    sdr->startup_delay+=1;
    fprintf(stderr, "   startup_delay=%i\n",sdr->startup_delay);
    return 0;
  }


  /* complex data conversion (also unsigned byte to signed) ; aka I/Q-samples */
  for (j=0; j<DAB_T_FRAME*2; j+=2) {
    sdr->real[j/2] = sdr->buffer[j] - 127; // = SCHAR_MAX
    sdr->imag[j/2] = sdr->buffer[j+1] - 127;
  }

  /* coarse time sync */
  /* performance bottleneck atm */
  sdr->coarse_timeshift = dab_coarse_time_sync(sdr->real,sdr->filt,sdr->force_timesync);
  // we are not in sync so -> next frame
  sdr->force_timesync=0;
  if (sdr->coarse_timeshift != 0) {
    //fprintf(stderr, "   coarse time shift %d\n", sdr->coarse_timeshift);
    return 0;
  }

  /* create complex frame */
  for (j=0; j<DAB_T_FRAME; j++) {
    sdr->dab_frame[j][0] = sdr->real[j];
    sdr->dab_frame[j][1] = sdr->imag[j];
  }

  /* fine time sync */
  sdr->fine_timeshift = dab_fine_time_sync(sdr->dab_frame);


  /* coarse_frequency shift */
  // fftw_plan fftw_plan_dft_1d(int n, fftw_complex *in, fftw_complex *out, int sign, unsigned flags);
  fftw_plan p;
  // FIXME why 505=DAB_T_GUARD+1??
  p = fftw_plan_dft_1d(DAB_T_CS, 
		       &sdr->dab_frame[DAB_T_NULL + DAB_T_GUARD/2], 
		       sdr->symbols[0], 
		       FFTW_FORWARD, 
		       FFTW_ESTIMATE);
  fftw_execute(p);
  fftw_destroy_plan(p);
  swap_blocks(DAB_T_CS, sdr->symbols[0]); // swap upper and lower blocks (of the PRS)

  int32_t coarse_freq_shift = dab_coarse_freq_sync_2(sdr->symbols[0]);
  if (abs(coarse_freq_shift) > 1) {
    // note: don't do the coarse shift immediately, it tends to irreversbly run away with the frequency
    // this is handled in demod_thread_fn
    sdr->coarse_freq_shift = coarse_freq_shift;
    sdr->force_timesync = 1;
    //fprintf(stderr, "   coarse freq shift, also force timeshift %d\n", sdr->coarse_freq_shift);
    // (even if req shift not applied immediately) return anyway to do the forced time sync
    return 0;
  }

  /* fine freq correction */
  sdr->fine_freq_shift = dab_fine_freq_corr(sdr->dab_frame);

  /* d-qpsk */
  // see https://en.wikipedia.org/wiki/Phase-shift_keying  (differential quad)
  // processing i=0 (again) is needed if usage (phase) of DAB_T_GUARD above
  // is different from here, because that changes relative phase too much already!
  for (i=0; i<DAB_SYMBOLS_IN_FRAME; i++) {
    p = fftw_plan_dft_1d(DAB_T_CS, 
			 &sdr->dab_frame[DAB_T_NULL + (DAB_T_SYM*i) + DAB_T_GUARD/2],
			 sdr->symbols[i], 
			 FFTW_FORWARD,
			 FFTW_ESTIMATE);
    fftw_execute(p);
    fftw_destroy_plan(p);
    swap_blocks(DAB_T_CS, sdr->symbols[i]);
  }
  for (j=1; j<DAB_SYMBOLS_IN_FRAME; j++) {
    for (i=0; i<DAB_T_CS; i++) {
      dqpsk_step(sdr->symbols[j][i], sdr->symbols[j-1][i], &sdr->symbols_d[j*DAB_T_CS + i]);
    }
  }
  
  uint8_t* dst = tf->fic_symbols_demapped[0];
  tf->has_fic = 1;  /* Always true for SDR input */

  int k,kk;
  for (j=1; j<DAB_SYMBOLS_IN_FRAME; j++) {
    if (j == 4) { dst = tf->msc_symbols_demapped[0]; }
    k = 0;
    for (i=0; i<DAB_T_CS; i++) {
      int empty_channels = DAB_T_CS/8;
      // 256..1792
      if (empty_channels <= i && i != DAB_T_CS/2 && i <= DAB_CARRIERS+empty_channels) {
        kk = rev_freq_deint_tab[k++];
        dst[kk]              = (sdr->symbols_d[j*DAB_T_CS+i][0] < 0);
        dst[kk+DAB_CARRIERS] = (sdr->symbols_d[j*DAB_T_CS+i][1] > 0);
      }
    }
    dst += 2*DAB_CARRIERS;  // 3072
  }
  
  //fprintf(stderr, "   ok\n");
  return 1;
}



void sdr_init(struct sdr_state_t *sdr)
{
  // circular buffer init
  cbInit(&(sdr->fifo), (DAB_T_FRAME*2*4)); // 4 frames; what is the 2? real/complex?
  // malloc of various buffers
  sdr->coarse_timeshift = 0;
  sdr->fine_timeshift=0;

  // malloc of various buffers
  sdr->symbols_d = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * DAB_T_CS * DAB_SYMBOLS_IN_FRAME);

  /* make sure to disable fault injection by default */
  sdr->p_e_prior_dep = 0.0f;
  sdr->p_e_prior_vitdec = 0.0f;
  sdr->p_e_after_vitdec = 0.0f;
}
