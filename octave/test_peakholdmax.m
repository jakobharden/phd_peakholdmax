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
  Fs = 4e3; # Hz, sampling rate
  Nc = 2; # number of cycles
  F = 1; # signal frequency
  spar = [1, F]; # A = 1 V, F = 1 Hz
  
  ## noise parameters
  Pv = 1; # V^2, noise power
  Nmc = 100; # number of Monte-Carlo simulation turns
  SNR = 20; # dB, signal-to-noise ratio
  
  ## generate noise standard
  N = floor(Nc * Fs / spar(2));
  vv_mat = tool_gen_noise(N, Nmc, Pv, rpath);
  
  ## generate exemplary signals
  [ds, fh] = test_peakholdmax_signals(Fs, Nc, spar, SNR, vv_mat);
  test_peakholdmax_save(ds, fh, rpath, 'sig');
  return;
  
##  ## find detection limits w.r.t. the signal frequency
##  ## given: peak hold length
##  ## unknown: signal frequency (variation variable)
##  [ds, fh] = test_peakholdmax_frqvar();
##  test_peakholdmax_save(ds, fh, rpath, 'frqvar');
##  
  ## find detection limits w.r.t. the signal-to-noise ratio (SNR)
  ## given: signal frequency
  ## given: peak hold length
  ## unknown: signal-to-noise ratio
  [ds, fh] = test_peakholdmax_snrvar(Fs, Nc, spar, vv_mat(:, 1:10));
  test_peakholdmax_save(ds, fh, rpath, 'snrvar');
  
endfunction


function [r_ds, r_fh] = test_peakholdmax_signals(p_fs, p_nc, p_sp, p_snr, p_vm);
  ## Plot exemplary signals
  ##
  ## p_fs  ... sampling rate, Hz, <uint>
  ## p_nc  ... number of cycles, <uint>
  ## p_sp  ... signal parameter array [A, F], [<dbl>, <dbl>]
  ## p_snr ... signal-to-noise ratio, dB, <dbl>
  ## p_vm  ... noise standard matrix, [[<dbl>]]
  ## r_ds  ... return: analysis result structure, <struct>
  ## r_fh  ... return: figure handle, [<uint>]
  
  ## generate signal, row vector
  [ss, nn, N] = tool_gen_signal(p_fs, p_nc, p_sp);
  
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


function [r_ds, r_fh] = test_peakholdmax_snrvar(p_fs, p_nc, p_sp, p_vm);
  
  r_ds = [];
  
  ## generate signal, row vector
  [ss, nn, N] = tool_gen_signal(p_fs, p_nc, p_sp);
  Ps = meansq(ss); # signal poser
  
  ## number of variations
  Nsnrvar = 20;
  Nhplvar = 50;
  
  ## number of Monte-Carlo simulation turns
  Nmc = size(p_vm, 2);
  
  ## SNR array, dB
  snr1 = 0; snr2 = 63;
  snr_arr = linspace(0, 63, Nsnrvar);
  
  ## hold peak counter limit array, number of samples
  hpc1 = 1; hpc2 = p_fs;
  hpc_arr = fix(linspace(hpc1, hpc2, Nhplvar));
  
  ## solutions
  nmax_exact = floor(p_fs / 4);
  
  ## detection window, used to validate correct detction results
  det_win = [floor(p_fs / 8), floor(3 * p_fs / 8)];
  
  det_est_xx = linspace(-1, 1, Nhplvar);
  det_est_xx1 = linspace(-0.47, 0.47, Nhplvar) + 0.5;
  det_est_yy = (1 - det_est_xx .^ 8) .* 53;
  
  ## loop over SNR array
  det_mat = ones(Nsnrvar, Nhplvar);
  nmax_mat = zeros(Nsnrvar, Nhplvar, Nmc);
  for k1 = 1 : numel(snr_arr)
    for k2 = 1 : numel(hpc_arr)
      ## Monte-Carlo simulation
      for k3 = 1 : Nmc
        ## noise array, column vector
        [vv, ~, ~, ~] = tool_scale_noise2snr(Ps, p_vm(:, k3), snr_arr(k1));
        ## signal in noise, column vector
        xx = ss(:) + vv;
        ## detect local maximum
        [v_max, n_max] = tool_peakholdmax(xx, 1, hpc_arr(k2), nmax_exact + hpc2 + 10);
        ## save detection state
        det_state = (n_max >= det_win(1)) && (n_max <= det_win(2));
        det_mat(k1, k2) = min([det_mat(k1, k2), det_state]);
        ## save detection results
        if (det_state)
          nmax_mat(k1, k2, k3) = n_max;
        else
          nmax_mat(k1, k2, k3) = NaN;
        endif
      endfor
    endfor
  endfor
  
  ## compute Monte-Carlo simulation statistics
  mc_med = median(nmax_mat, 3);
  
  ## compute Monte-Carlo simulation statistics
  mc_std = std(nmax_mat, 0, 3);
  
  ## plot state
  r_fh(1) = figure('name', 'state', 'position', [100, 100, 800, 0.62*800]);
  ah = axes(r_fh(1), 'tickdir', 'out');
  hold(ah, 'on');
  imagesc(hpc_arr / p_fs, flip(snr_arr), det_mat);
  plot(det_est_xx1, det_est_yy, 'linewidth', 2, 'color', [1, 0, 0]);
  hold(ah, 'off');
  colorbar();
  title(ah, sprintf('Detection state'));
  xlabel(ah, 'HPC / N [1]');
  ylabel(ah, 'SNR [dB]');
  
  ## plot median
  r_fh(2) = figure('name', 'median', 'position', [100, 100, 800, 0.62*800]);
  ah = axes(r_fh(2), 'tickdir', 'out');
  hold(ah, 'on');
  imagesc(hpc_arr / p_fs, flip(snr_arr), mc_med);
  plot(det_est_xx1, det_est_yy, 'linewidth', 2, 'color', [1, 0, 0]);
  hold(ah, 'off');
  colorbar();
  title(ah, sprintf('Local maximum - median'));
  xlabel(ah, 'HPC / N [1]');
  ylabel(ah, 'SNR [dB]');
  
  ## plot empirical variance
  r_fh(3) = figure('name', 'empstd', 'position', [100, 100, 800, 0.62*800]);
  ah = axes(r_fh(3), 'tickdir', 'out');
  hold(ah, 'on');
  imagesc(hpc_arr / p_fs, flip(snr_arr), mc_std);
  plot(det_est_xx1, det_est_yy, 'linewidth', 2, 'color', [1, 0, 0]);
  hold(ah, 'off');
  colorbar();
  title(ah, sprintf('Local maximum - empirical standard deviation'));
  xlabel(ah, 'HPC / N [1]');
  ylabel(ah, 'SNR [dB]');
  
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
  save('-binary', fullfile(p_rp, 'oct', sprintf('%s.oct', p_pfx)), 'ads');
  
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
