## Generate Gaussian white noise matrix and store it in a file (standard noise that allows to reproduce results)
##
## Usage: [r_vv] = tool_gen_noise(p_nsmp, p_nsig, p_vpow, p_ofp)
##
## p_nsmp ... number of samples, rows of noise matrix, <uint>
## p_nsig ... number of signals, columns of noise matrix, <uint>
## p_vpow ... noise power, optional, default = 1, <dbl>
## p_ofp  ... output file path, optional, default = "./results", <str>
## r_vv   ... return: standard noise matrix, [Nsmp x Nsig], [[<dbl>]]
##
## Note: Existing file are not overwritten by this function. If a new noise standard is required, remove the existing noise standard file.
##
## see also: randn
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
function [r_vv] = tool_gen_noise(p_nsmp, p_nsig, p_vpow, p_ofp)
  
  ## check arguments
  if (nargin < 4)
    p_ofp = './results';
  endif
  if (nargin < 3)
    p_vpow = [];
  endif
  if (nargin < 2)
    help tool_gen_noise;
    error('Less arguments given!');
  endif
  if isempty(p_vpow)
    p_vpow = 1;
  endif
  
  ## standard noise file path
  nfp = fullfile(p_ofp, sprintf('noise_Nsmp%d_Nmc%d.oct', p_nsmp, p_nsig));
  
  ## check if file exists
  if (exist(nfp) == 3)
    r_vv = load(nfp, 'r_vv').r_vv;
    printf('tool_gen_noise: Loaded existing standard noise! File = %s\n', nfp);
  else
    ## generate noise
    if (p_vpow == 1)
      r_vv = randn(p_nsmp, p_nsig);
    else
      r_vv = randn(p_nsmp, p_nsig) .* sqrt(p_vpow);
    endif
    ## save noise standard
    save('-binary', nfp, 'r_vv');
    ## print message
    printf('tool_gen_noise: Saved new noise standard to %s. Variable name = r_vv\n', nfp);
  endif
  
endfunction
