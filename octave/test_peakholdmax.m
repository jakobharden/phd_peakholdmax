## Test the peak-hold-max algorithm
##
## Usage: test_peakholdmax()
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
function test_peakholdmax()
  
  ## results path
  rpath = './results/test_peakholdmax';
  
  ## signal parameters
  Fs = 1e4; # Hz, sampling rate
  Nc = 2; # number of cycles
  F = 1; # signal frequency
  spar = [1, F]; # A = 1 V, F = 1 Hz
  
  ## noise parameters
  Pv = 1; # V^2, noise power
  Nsig = 100; # number of signals
  Nv = Fs; # prefix section length, noise floor
  SNR = 20; # dB, signal-to-noise ratio
  
  ## generate noise standard
  N = floor(Nc * Fs / spar(2) + Nv);
  vv_mat = tool_gen_noise(N, Nsig, Pv, rpath);
  
  ## generate exemplary signals
  [ds, fh] = test_peakholdmax_signals(Fs, Nc, spar, Nv, SNR, vv_mat);
  test_peakholdmax_save(ds, fh, rpath, 'sig');
  
##  ## find detection limits w.r.t. the signal frequency
##  ## given: peak hold length
##  ## unknown: signal frequency (variation variable)
##  [ds, fh] = test_peakholdmax_frqvar();
##  test_peakholdmax_save(ds, fh, rpath, 'frqvar');
##  
##  ## find detection limits w.r.t. the signal-to-noise ratio (SNR)
##  ## given: signal frequency
##  ## given: peak hold length
##  ## unknown: signal-to-noise ratio
##  [ds, fh] = test_peakholdmax_snrvar();
##  test_peakholdmax_save(ds, fh, rpath, 'snrvar');
  
endfunction


function [r_ds, r_fh] = test_peakholdmax_signals(p_fs, p_nc, p_sp, p_nv, p_snr, p_vm);
  ## Plot exemplary signals
  ##
  ## p_fs  ... sampling rate, Hz, <uint>
  ## p_nc  ... number of cycles, <uint>
  ## p_nv  ... prefix section length, noise floor, <uint>
  ## p_sp  ... signal parameter array [A, F], [<dbl>, <dbl>]
  ## p_snr ... signal-to-noise ratio, dB, <dbl>
  ## p_vm  ... noise standard matrix, [[<dbl>]]
  ## r_ds  ... return: analysis result structure, <struct>
  ## r_fh  ... return: figure handle, [<uint>]
  
  ## generate signal, row vector
  [ss, nn, N] = tool_gen_signal(p_fs, p_nc, p_sp);
  ss = [zeros(1, p_nv), ss];
  N = N + p_nv;
  nn = 1 : N;
  
  ## get noise, column vector
  [vv, ~, ~, ~] = tool_scale_noise2snr(ss, p_vm(:, 1), p_snr);
  
  ## generate signal corrupted by noise, column vector
  xx = ss(:) + vv;
  
  ## plot signal
  r_fh = figure('name', 'Signal', 'position', [100, 100, 800, 0.62*800]);
  ah = axes(r_fh, 'tickdir', 'out');
  hold(ah, 'on');
  plot(ah, nn, xx, ';x[n];', 'color', ones(1, 3) * 0.5, 'linewidth', 0.5);
  plot(ah, nn, ss, ';s[n];', 'color', [1, 0, 0], 'linewidth', 2);
  hold(ah, 'off');
  title(ah, sprintf('Test signal\nF_s = %d Hz, A = %d V, F = %d Hz, SNR = %d dB', p_fs, p_sp(1), p_sp(2), p_snr));
  xlabel(ah, 'Signal index [1]');
  ylabel(ah, 'Amplitude [V]');
  
endfunction


function [r_ds, r_fh] = test_peakholdmax_frqvar();
  
  
endfunction


function [r_ds, r_fh] = test_peakholdmax_snrvar();
  
  
endfunction


function [r_ds, r_fh] = test_peakholdmax_save(p_ds, p_fh, p_rp, p_pfx);
  ## Save results
  ##
  ## p_ds  ... analysis result data structure, <struct>
  ## p_fh  ... figure handle(s), [<uint>]
  ## p_rp  ... result directory path, <str>
  ## p_pfx ... file name prefix, <str>
  
  ## export binary data file
  ads = p_ds;
  save('-binary', fullfile(p_rp, 'oct', sprintf('%s.oct', fnpfx)), 'ads');
  
  ## export figures
  for j = 1 : length(p_fh)
    if (length(p_fh) == 1)
      idxstr = '';
    else
      idxstr = sprintf('_%d', j);
    endif
    ## save figure
    hgsave(p_fh(j), fullfile(p_rp, 'fig', sprintf('%s%s.ofig', p_pfx, idxstr)));
    ## save bitmap
    saveas(p_fh(j), fullfile(p_rp, 'png', sprintf('%s%s.png', p_pfx, idxstr)), 'png');
    ## close figure
    close(p_fh(j));
  endfor
  
endfunction
