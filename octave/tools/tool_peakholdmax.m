## Detect the first local maximum in a signal amplitude array (peak-hold detection)
##
## p_xx  ... signal amplitude array, [<dbl>]
## p_n0  ... detection window start index, <uint>
## p_cl  ... hold-peak counter limit, <uint>
## p_dl  ... detection window length, optional, <uint>
## r_max ... return: local maximum value, <dbl>
## r_idx ... return: local maximum index, <uint>
##
## Requirements:
##   1) discrete sampled low-pass signal in noise, equally spaced samples
##   2) the detection window must contain at least one maximum
##   3) the amplitude of the first local maximum must be greater than 20 percent of the global maximum amplitude
##   4) the hold-peak counter limit p_cl must be lower than the smallest distance between the local maxima of the signal
##
## Algorithm:
##   step 1: determine global maximum value and sample index (initialize return values)
##   step 2: determine maximum detection window length dl_max
##   step 3: check detection window length, condition: p_dl == []; true: p_dl = dl_max; false: p_dl = min([p_dl, dl_max])
##   step 4: replace signal amplitudes lower than 20% of global maximum amplitude with zero (noise floor removal)
##   step 5: initialize loop variables, a = temporary amplitude maximum, c = hold-peak counter
##   step 6: loop over signal sample indices n within the detection window
##     step 6.1: evaluate update-condition, condition: x[n] >= a; true: a= x[n], c = 0; false: c += 1
##     step 6.2: evaluate stop-condition, condition: c == p_cl; true: r_max = a, r_idx = p_cl - n, break looping; false: next loop
##
#######################################################################################################################
## LICENCE
##
##    Copyright (C) 2025 Jakob Harden (jakob.harden@tugraz.at, Graz University of Technology, Graz, Austria)
##    This file is part of the PhD thesis of Jakob Harden.
##
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU Affero General Public License as
##    published by the Free Software Foundation, either version 3 of the
##    License, or (at your option) any later version.
##
##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU Affero General Public License for more details.
##
##    You should have received a copy of the GNU Affero General Public License
##    along with this program.  If not, see <https://www.gnu.org/licenses/>.
##
#######################################################################################################################
##
function [r_max, r_idx] = tool_peakholdmax(p_xx, p_n0, p_cl, p_dl)
  
  ## check arguments
  if (nargin < 4)
    p_dl = [];
  endif
  if (nargin < 3)
    help tool_peakholdmax;
    error('Less arguments given!');
  endif
  
  ## step 1: determine global maximum value and sample index (init return values)
  [r_max, k_globmax] = max(p_xx(p_n0 : end)); # global maximum value and local sample index
  r_idx = k_globmax + p_n0 - 1; # map local to global sample index
  
  ## step 2: determine maximum detection window length
  dl_max = numel(p_xx) - p_n0 + 1;
  
  ## step 3: check detection window length
  if isempty(p_dl)
    ## step 3.1: use maximum window length
    p_dl = dl_max;
  else
    ## step 3.2: limit detection window length to maximum window length
    p_dl = min([p_dl, dl_max]);
  endif
  
  ## step 4: replace signal amplitudes lower than 20% of global maximum amplitude with zero
  ## this is a measure to reduce the likelyhood of noise-peak detection (noise floor removal)
  rm_idx = find(abs(p_xx) < abs(r_max) * 0.2);
  p_xx(rm_idx) = 0.0;
  
  ## step 5: initialize loop variables
  a = min(p_xx(p_n0 : end)); # temporary maximum amplitude, use global minimum
  c = 0; # hold-peak counter
  
  ## step 6: loop over signal samples within detection window
  for n = p_n0 : (p_n0 + p_dl - 1)
    ## step 6.1: evaluate update-condition, x[n] >= a
    if (p_xx(n) >= a)
      ## update running maximum amplitude and hold-peak counter
      a = p_xx(n);
      c = 0;
    else
      ## increment hold-peak counter only
      c += 1;
    endif
    ## step 6.2: evaluate stop-condition, c == p_cl
    if (c == p_cl)
      ## reached counter limit, update return values
      r_max = a;
      r_idx = n - p_cl;
      break;
    endif
  endfor
  
endfunction
