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

#include <stdio.h>
#include <stdint.h>
#include <math.h>
#include <fftw3.h>
#include <stdlib.h>
#include "dab_constants.h"

double cpx_abs(fftw_complex sample);

int time_sync_still_good_2(fftw_complex * frame);
int time_sync_full(fftw_complex * frame);
int dab_fine_time_sync(fftw_complex * frame);

int dab_coarse_freq_sync_4(int symcount, fftw_complex (*symbols)[DAB_T_CS]);
double dab_fine_freq_corr(fftw_complex * dab_frame);
double ffs_from_phase(fftw_complex * symbols_d);


