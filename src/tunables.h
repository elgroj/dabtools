#ifndef _TUNABLES_H
#define _TUNABLES_H
// magic numbers for performance/reliability tuning


// full usage of viterbi: SCALEs=128, GAIN=4.0
//  -> takes 75% of the time (and doesn't seem to be so much better???
// originally behaviour: both SCALE=1; GAIN high (100 or so)
//  -> this is actually only faster if it fails a lot!!!

// orig value: 1.0
#define VITERBI_NOISE 0.7

// int; max for char is 128; min is 1; 
#define VITERBI_SCALE 128
// scale where to clamp the values; 128 uses the full byte
#define VITERBI_CLAMP_SCALE 128
// obtained by experimentation: 4.0; 
// extra gain (besides VITERBI_SCALE) for converting normalized dqpsk values into viterbi symbols
#define SCALING_GAIN 4.0

// magic threshold was 50.0
// fine frequency shifts below this threshold in Hz will be ignored
// (there is always a risk of loosing lock after a frequency change)
#define FFS_THRESHOLD 35.0
// number of frames to ignore frequency shifts, so things can stabilize
#define FS_BLOCKING 4

// number of symbols to check for coarse frequency shift
#define CFS_SYM_COUNT 15

// sampling skip used for coarse time sync 
// number of samples 
//   in check: DAB_T_NULL/FILT_SAMPLING
//   at sync: (DAB_T_FRAME-DAB_T_NULL)/FILT_SAMPLING + 1
#define FILT_SAMPLING 10

// threshold factor by which the null symbol should be more quiet than the data frame
#define REALATIVE_ENERGY_THRESHOLD 0.2

// similar thing for new time sync: first pass will be done with this as a stride
#define TIME_SYNC_INIT_STEP 5
#define FINE_TIME_SHIFT_MAX 1000

#endif
