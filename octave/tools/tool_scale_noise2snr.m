## Scale noise amplitudes to given signal-to-noise ratio (SNR)
##
## Usage: [r_vv, r_sf, r_ps, r_pv] = tool_scale_noise2snr(p_ss, p_vv, p_snr)
##
## p_ss  ... signal amplitude array or signal power, [<dbl>] or <dbl>
## p_vv  ... noise amplitude array, [<dbl>]
## p_snr ... requested signal-to-noise ratio (SNR), dB, <dbl>
## r_vv  ... return: scaled noise amplitude array, [<dbl>]
## r_sf  ... return: noise amplitude scaling factor, <dbl>
## r_ps  ... return: signal power, <dbl>
## r_pv  ... return: power of scaled noise amplitude array, <dbl>
##
## see also: meansq
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
function [r_vv, r_sf, r_ps, r_pv] = tool_scale_noise2snr(p_ss, p_vv, p_snr)
  
  ## check arguments
  if (nargin < 3)
    help tool_scale_noise2snr;
    error('Less arguments given!');
  endif
  
  ## init return values
  r_vv = [];
  r_sf = 0;
  r_ps = 0;
  r_pv = 0;
  
  ## signal power
  if (length(p_ss) < 2)
    ## given: signal power
    r_ps = p_ss;
  else
    ## given: signal amplitude array
    r_ps = meansq(p_ss);
  endif
  
  ## noise power, unscaled
  pv = meansq(p_vv);
  
  ## check SNR
  if isinf(p_snr)
    r_vv = zeros(size(p_vv));
    return;
  endif
  
  ## noise power w.r.t. SNR
  ## SNR = 10 * log10(ps / pv) => pv_snr = ps / 10^(SNR / 10)
  p_snr = abs(p_snr);
  pv_snr = r_ps / (10^(p_snr / 10));
  r_pv = pv_snr;
  
  ## noise amplitude scaling factor
  r_sf = sqrt(pv_snr / pv);
  
  ## scale noise amplitudes
  r_vv = p_vv * r_sf;
  
endfunction
