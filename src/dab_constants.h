#ifndef _DAB_CONSTANTS_H
#define _DAB_CONSTANTS_H

// L 76       number of OFDM symbols per transmission frame (the Null symbol being excluded);
#define DAB_SYMBOLS_IN_FRAME 76

// K 1536                   number of transmitted carriers
#define DAB_CARRIERS 1536

// T_F  196608 T     96 ms  transmission frame duration
#define DAB_T_FRAME 196608

// T_NULL 2656 T ~1,297 ms  Null symbol duration
#define DAB_T_NULL 2656

// T_s    2552 T ~1,246 ms  duration of OFDM symbols of indices l = 1..L
#define DAB_T_SYM 2552

// T_u    2048 T  1     ms  inverse of the carrier spacing
#define DAB_T_CS 2048

// ∆       504 T   ~246 μs  duration of the time interval called guard interval
#define DAB_T_GUARD 504

//  checked in sdr_demod:
//  assert(DAB_T_FRAME == DAB_SYMBOLS_IN_FRAME * DAB_T_SYM + DAB_T_NULL);
//  assert(DAB_T_SYM == DAB_T_CS + DAB_T_GUARD);

#endif
