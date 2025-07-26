## Test the local maximum detection algorithm
##
## Usage: test_peakholdmax(p_rm)
##
## p_rm ... run mode, optional, default = 'compileplot', <str>
##
## Run modes:
##   p_rm = 'compile': compile the analysis results
##   p_rm = 'plot': plot already compiled analysis results
##   p_rm = 'compileplot': compile and plot the analysis results, default
##
## Note: All analysis results and plots are stored in the "./results/test_peakholdmax/" directory.
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
function test_peakholdmax(p_rm = 'compileplot')
  
  ## results path
  rpath = './results/test_peakholdmax';
  
  ## signal parameters
  Fs = 4e3; # Hz, sampling frequency
  Nc = 2; # number of cycles
  A = 1; # amplitude
  F = 1; # signal frequency
  SNR = 20; # dB, signal-to-noise ratio
  decay_arr = [-0.50, -0.25, 0.00, 0.25, 0.50]; # damping/amplification factors
  
  ## variation parameters
  Nmc = 500; # number of Monte-Carlo simulation turns
  Nclimv = 21; # counter limit variation number (21 --> 5 % steps)
  Nsnrv = 22; # SNR variation number (22 --> 3 dB steps)
  
  ## plot parameters
  cm_gradient = 'broc'; # colour map, see also ./scicol_crameri/colormap_crameri.m
  
  ## generate noise standard matrix (this allows to reproduce the analysis results)
  Nvv = floor(Nc * Fs / F);
  vv_mat = tool_gen_noise(Nvv, Nmc, Pv = 1, rpath);
  
  ## generate test signals
  for k = 1 : numel(decay_arr)
    spar = [A, F, decay_arr(k)]; # amplitude, frequency, exponential decay factor
    pfx = sprintf('sig%d', k);
    switch (lower(p_rm))
      case 'compile'
        ds = test_peakholdmax_signals_compile(Fs, Nc, spar, SNR, vv_mat);
        test_peakholdmax_save_analysis(ds, rpath, pfx);
      case 'plot'
        ds = test_peakholdmax_load_analysis(rpath, pfx);
        fh = test_peakholdmax_signals_plot(ds);
        test_peakholdmax_save_figure(fh, rpath, pfx);
      otherwise
        ds = test_peakholdmax_signals_compile(Fs, Nc, spar, SNR, vv_mat);
        test_peakholdmax_save_analysis(ds, rpath, pfx);
        fh = test_peakholdmax_signals_plot(ds);
        test_peakholdmax_save_figure(fh, rpath, pfx);
    endswitch
  endfor

  ## determine detection limits w.r.t. the signal-to-noise ratio (SNR) and the exponential decay factor
  for k = 1 : numel(decay_arr)
    spar = [A, F, decay_arr(k)]; # amplitude, frequency, decay factor
    pfx = sprintf('snrvar%d', k);
    switch (lower(p_rm))
      case 'compile'
        ds = test_peakholdmax_snrvar_compile(Fs, spar, vv_mat, Nclimv, Nsnrv);
        test_peakholdmax_save_analysis(ds, rpath, pfx)
      case 'plot'
        ds = test_peakholdmax_load_analysis(rpath, pfx);
        fh = test_peakholdmax_snrvar_plot(ds, cm_gradient);
        test_peakholdmax_save_figure(fh, rpath, pfx);
      otherwise
        ds = test_peakholdmax_snrvar_compile(Fs, spar, vv_mat, Nclimv, Nsnrv);
        test_peakholdmax_save_analysis(ds, rpath, pfx)
        fh = test_peakholdmax_snrvar_plot(ds, cm_gradient);
        test_peakholdmax_save_figure(fh, rpath, pfx);
    endswitch
  endfor
  
endfunction


function [r_ds] = test_peakholdmax_signals_compile(p_fs, p_nc, p_sp, p_snr, p_vm)
  ## Compile data for exemplary signals
  ##
  ## p_fs  ... sampling frequency, Hz, <uint>
  ## p_nc  ... number of cycles, <uint>
  ## p_sp  ... signal parameter array [A, F, beta], [<dbl>, <dbl>, <dbl>]
  ## p_snr ... signal-to-noise ratio, dB, <dbl>
  ## p_vm  ... noise standard matrix, [[<dbl>]]
  ## r_ds  ... return: analysis result structure, <struct_signal>
  
  ## generate signal, row vector
  [ss, nn, N] = tool_gen_signal(p_fs, p_nc, p_sp);
  
  ## determine exact maximum point
  [v_max, n_max] = max(ss(1 : floor(N / (2 * p_nc))));
  
  ## window of accepted solutions, window length is 25 % of the wave length
  dl = round(p_fs / (8 * p_sp(2)));
  wd = [n_max - dl, n_max + dl - 1];
  
  ## get noise, column vector
  [vv, ~, ~, ~] = tool_scale_noise2snr(ss, p_vm(:, 1), p_snr);
  
  ## generate signal corrupted by noise, column vector
  xx = ss(:) + vv;
  
  ## detected maximum points
  [v_maxd1, n_maxd1] = max(xx(1 : 4 * dl)); # accepted detection
  [v_maxd2, n_maxd2] = max(xx(8 * dl : 12 * dl));   n_maxd2 = n_maxd2 + 8 * dl - 1; # false detection
  
  ## create data structure
  r_ds.obj = 'struct_signal';
  r_ds.ver = uint16([1, 0]);
  r_ds.des = 'Signal data structure';
  r_ds.author_surename = 'Harden';
  r_ds.author_givenname = 'Jakob';
  r_ds.author_affiliation = 'Graz University of Technology, Graz, Austria';
  r_ds.author_orcid = 'https://orcid.org/0000-0002-5752-1785';
  r_ds.author_email1 = 'jakob.harden@tugraz.at';
  r_ds.author_email2 = 'jakob.harden@student.tugraz.at';
  r_ds.author_email3 = 'office@jakobharden.at';
  r_ds.author_website = 'https://jakobharden.at/wordpress/';
  r_ds.licence = 'This file is licensed under a Creative Commons Attribution 4.0 International license';
  r_ds.p_fs = p_fs;
  r_ds.p_nc = p_nc;
  r_ds.p_sp = p_sp;
  r_ds.p_snr = p_snr;
  r_ds.ss = ss;
  r_ds.nn = nn;
  r_ds.N = N;
  r_ds.v_max = v_max;
  r_ds.n_max = n_max;
  r_ds.dl = dl;
  r_ds.wd = wd;
  r_ds.vv = vv;
  r_ds.xx = xx;
  r_ds.v_maxd1 = v_maxd1;
  r_ds.n_maxd1 = n_maxd1;
  r_ds.v_maxd2 = v_maxd2;
  r_ds.n_maxd2 = n_maxd2;
  
endfunction


function [r_ds] = test_peakholdmax_snrvar_compile(p_fs, p_sp, p_vm, p_ncv, p_nsv)
  ## Compile SNR variation, Monte-Carlo simulation
  ##
  ## p_fs  ... sampling frequency, Hz, <dbl>
  ## p_sp  ... signal parameter array [A, F, beta], [<dbl>, <dbl>, <dbl>]
  ## p_vm  ... noise standard matrix, [[<dbl>]]
  ## p_ncv ... number of c_lim variations, <uint>
  ## p_nsv ... number of SNR variations, <uint>
  ## r_ds  ... return: analysis result data structure, <struct_analysis>
  
  ## generate signal, row vector
  [ss, nn, N] = tool_gen_signal(p_fs, Ncy = 2, p_sp);
  Ps = meansq(ss); # signal power
  N1 = p_fs; # number of samples of one full wave
  
  ## exact maximum point
  [v_max_ex, n_max_ex] = max(ss(1 : floor(N / (2 * Ncy))));
  
  ## number of variations
  Nsnrvar = p_nsv; # SNR
  Nclimvar = p_ncv; # counter limit
  Ntotvar = Nsnrvar * Nclimvar; # total number of parameter variations
  
  ## number of Monte-Carlo simulation turns
  Nmc = size(p_vm, 2);
  
  ## SNR array, dB
  snr_arr = linspace(0, 63, Nsnrvar);
  d_snr = snr_arr(2) - snr_arr(1);
  
  ## hold peak counter limit array, number of samples
  clim_arr = round(linspace(1, N1, Nclimvar));
  d_clim = clim_arr(2) - clim_arr(1);
  
  ## window of accepted solutions, window length is 25 % of the wave length
  dl = round(p_fs / (8 * p_sp(2)));
  det_win = [n_max_ex - dl, n_max_ex + dl - 1];
  
  ## initialise analysis results
  det_mat = ones(Nsnrvar, Nclimvar);
  ed_mean = zeros(Nsnrvar, Nclimvar);
  ed_std = zeros(Nsnrvar, Nclimvar);
  er_mean = zeros(Nsnrvar, Nclimvar);
  er_std = zeros(Nsnrvar, Nclimvar);
  nmax_arr = zeros(Nmc, 1);
  
  ## loop over SNR array
  disp('Running Monte-Carlo simulation');
  for k1 = 1 : Nsnrvar
    ## loop over counter limit array
    for k2 = 1 : Nclimvar
      printf('  Variation %d of %d\n', (k1 - 1) * Nclimvar + k2, Ntotvar)
      ## Monte-Carlo simulation
      for k3 = 1 : Nmc
        ## noise array, column vector
        [vv, ~, ~, ~] = tool_scale_noise2snr(Ps, p_vm(:, k3), snr_arr(k1));
        ## signal in noise, column vector
        xx = ss(:) + vv;
        ## detect local maximum
        [v_max, n_max] = tool_det_localmax(xx, 1, clim_arr(k2), N);
        ## save detection state
        det_state = (n_max >= det_win(1)) && (n_max <= det_win(2));
        det_mat(k1, k2) = min([det_mat(k1, k2), det_state]);
        ## save detection results
        if (det_state)
          nmax_arr(k3) = n_max;
        else
          nmax_arr(k3) = NaN;
        endif
      endfor
      ## compile error and error statistics
      ed = nmax_arr - n_max_ex; # error, number of samples
      er = ed / N1 * 100; # error, relative to signal frequency, percent
      ed_mean(k1, k2) = mean(ed);
      ed_std(k1, k2) = std(ed);
      er_mean(k1, k2) = mean(er);
      er_std(k1, k2) = std(er);
    endfor
  endfor
  
  ## create data structure
  r_ds.obj = 'struct_analysis';
  r_ds.ver = uint16([1, 0]);
  r_ds.des = 'Analysis data structure';
  r_ds.author_surename = 'Harden';
  r_ds.author_givenname = 'Jakob';
  r_ds.author_affiliation = 'Graz University of Technology, Graz, Austria';
  r_ds.author_orcid = 'https://orcid.org/0000-0002-5752-1785';
  r_ds.author_email1 = 'jakob.harden@tugraz.at';
  r_ds.author_email2 = 'jakob.harden@student.tugraz.at';
  r_ds.author_email3 = 'office@jakobharden.at';
  r_ds.author_website = 'https://jakobharden.at/wordpress/';
  r_ds.licence = 'This file is licensed under a Creative Commons Attribution 4.0 International license';
  r_ds.ss = ss;
  r_ds.nn = nn;
  r_ds.N = N;
  r_ds.Ncy = Ncy;
  r_ds.Ps = Ps;
  r_ds.N1 = N1;
  r_ds.v_max_ex = v_max_ex;
  r_ds.n_max_ex = n_max_ex;
  r_ds.Nsnrvar = Nsnrvar;
  r_ds.Nclimvar = Nclimvar;
  r_ds.Ntotvar = Ntotvar;
  r_ds.Nmc = Nmc;
  r_ds.snr_arr = snr_arr;
  r_ds.d_snr = d_snr;
  r_ds.clim_arr = clim_arr;
  r_ds.d_clim = d_clim;
  r_ds.dl = dl;
  r_ds.det_win = det_win;
  r_ds.det_mat = det_mat;
  r_ds.ed_mean = ed_mean;
  r_ds.ed_std = ed_std;
  r_ds.er_mean = er_mean;
  r_ds.er_std = er_std;
  
endfunction


function test_peakholdmax_save_analysis(p_ds, p_rp, p_pfx)
  ## Save analysis results
  ##
  ## p_ds  ... analysis result data structure, <struct>
  ## p_rp  ... result directory path, <str>
  ## p_pfx ... file name prefix, <str>
  
  ## export binary data file
  ads = p_ds;
  save('-binary', fullfile(p_rp, 'oct', sprintf('%s.oct', p_pfx)), 'ads');
  
endfunction


function [r_ds] = test_peakholdmax_load_analysis(p_rp, p_pfx)
  ## Load analysis results
  ##
  ## p_rp  ... result directory path, <str>
  ## p_pfx ... file name prefix, <str>
  ## r_ds  ... return: analysis result data structure, <struct>
  
  ## load binary data file
  r_ds = load(fullfile(p_rp, 'oct', sprintf('%s.oct', p_pfx)), 'ads').ads;
  
endfunction


function test_peakholdmax_save_figure(p_fh, p_rp, p_pfx)
  ## Save figures
  ##
  ## p_fh  ... figure handle(s), [<uint>]
  ## p_rp  ... result directory path, <str>
  ## p_pfx ... file name prefix, <str>
  
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


function [r_fh] = test_peakholdmax_signals_plot(p_ds)
  ## Plot exemplary signals
  ##
  ## p_ds  ... signal data structure, <struct_signal>
  ## r_fh  ... return: figure handle, <uint>
  
  ## plot signal
  y_max = max(p_ds.xx);
  y_min = min(p_ds.xx);
  r_fh = figure('name', 'Signal', 'position', [100, 100, 800, 0.62*800]);
  annotation(r_fh, 'textbox', [0.005, 0.005, 1, 0.1], 'fontsize', 8, ...
    'string', 'CC BY-4.0 (2025) Jakob Harden (Graz University of Technology)', 'linestyle', 'none');
  ah = axes(r_fh, 'tickdir', 'out');
  xlim(ah, [1 - p_ds.dl/10, p_ds.nn(end) + p_ds.dl/10]);
  hold(ah, 'on');
  fill(ah, [p_ds.wd, flip(p_ds.wd)], [y_min, y_min, y_max, y_max], [1.0000, 0.7137, 0.4235], 'handlevisibility', 'off', 'linestyle', 'none'); # window, accepted solution
  plot(ah, [p_ds.nn(1), p_ds.nn(end)], [0, 0], 'color', [1, 1, 1] * 0.25, 'handlevisibility', 'off'); # x axis
  plot(ah, p_ds.nn, p_ds.xx, ';x[n];', 'color', ones(1, 3) * 0.5, 'linewidth', 0.5); # noisy signal
  plot(ah, p_ds.nn, p_ds.ss, ';s[n];', 'color', [1, 0, 0], 'linewidth', 2); # clean signal
  plot(ah, 1, 0, 'ok', 'linewidth', 2, 'handlevisibility', 'off'); # start point
  plot(ah, p_ds.n_max, p_ds.v_max, 'ok', 'linewidth', 2, 'handlevisibility', 'off'); # exact maximum point
  plot(ah, [1, 1] * p_ds.n_max, [0, p_ds.v_max], ':k', 'linewidth', 1, 'handlevisibility', 'off'); # exact maximum point, line
  plot(ah, p_ds.n_maxd1, p_ds.v_maxd1, 'og', 'linewidth', 2, 'handlevisibility', 'off'); # detected maximum point, accepted solution
  if (p_ds.v_maxd2 >= p_ds.v_maxd1)
    plot(ah, p_ds.n_maxd2, p_ds.v_maxd2, 'or', 'linewidth', 2, 'handlevisibility', 'off'); # detected maximum point, false detection
  endif
  plot(ah, [p_ds.n_max, p_ds.n_max * 2.5], [1, 1] * p_ds.v_maxd1, '->k', 'handlevisibility', 'off'); # arrow, c_lim
  plot(ah, [p_ds.n_max * 2.5, p_ds.nn(end)], [1, 1] * p_ds.v_maxd1, '--k', 'handlevisibility', 'off'); # arrow extension, c_lim
  text(ah, 1, -p_ds.p_sp(1) * 0.05, ' n_i', 'fontsize', 12, 'verticalalignment', 'top');
  text(ah, p_ds.n_max, -p_ds.p_sp(1) * 0.05, ' n_{max,exact}', 'fontsize', 12, 'verticalalignment', 'top');
  text(ah, mean([p_ds.n_max, p_ds.n_max * 2.5]), p_ds.v_maxd1 * 0.98, 'c_{lim}', 'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'center');
  text(ah, mean(p_ds.wd), y_min, sprintf('window of\naccepted\nsolutions'), 'fontsize', 10, 'verticalalignment', 'bottom', 'horizontalalignment', 'center');
  text(ah, p_ds.n_maxd1, p_ds.v_maxd1 * 1.05, '  max(x), accepted', 'fontsize', 10, 'verticalalignment', 'bottom', 'horizontalalignment', 'center');
  if (p_ds.v_maxd2 >= p_ds.v_maxd1)
    text(ah, p_ds.n_maxd2, p_ds.v_maxd2 * 1.05, '  max(x), false', 'fontsize', 10, 'verticalalignment', 'bottom', 'horizontalalignment', 'center');
  endif
  hold(ah, 'off');
  title(ah, sprintf('Test signal\nF_s = %d Hz, A = %d V, F = %d Hz, E = %.2f, SNR = %d dB', p_ds.p_fs, p_ds.p_sp(1), p_ds.p_sp(2), p_ds.p_sp(3), p_ds.p_snr));
  xlabel(ah, 'Signal index n [1]');
  ylabel(ah, 'Amplitude A [V]');
  legend(ah, 'boxoff');
  
endfunction


function [r_fh] = test_peakholdmax_snrvar_plot(p_ds, p_cm)
  ## Plot SNR variation results
  ##
  ## p_ds ... analysis result data structure, <struct_analysis>
  ## p_cm ... color map name, see also ./scicol_crameri/colormap_crameri.m, <str>
  ## r_fh ... return: figure handle, <uint>
  
  ## select scientific color map
  cmap = colormap_crameri(p_cm);
  ## interpolate colour maps for error representation (zero = white colour)
  cmap2 = test_peakholdmax_intpcolmap(p_ds.ed_mean, cmap);
  cmap3 = test_peakholdmax_intpcolmap(p_ds.ed_std, cmap);
  cmap4 = test_peakholdmax_intpcolmap(p_ds.er_mean, cmap);
  cmap5 = test_peakholdmax_intpcolmap(p_ds.er_std, cmap);
  
  ## upper zone, SNR_lim1 <= SNR <= SNR_max
  SNR_lim1 = 24; # dB
  [clim_lo1, clim_hi1, rr1, cc1] = tool_peakholdmax_clim_wrt_snr(p_ds.det_mat, [SNR_lim1, p_ds.snr_arr(end)], p_ds.snr_arr, p_ds.clim_arr, p_ds.N1);
  
  ## lower zone, SNR_lim2 <= SNR <= SNR_lim1
  SNR_lim2 = 12; # dB
  [clim_lo2, clim_hi2, rr2, cc2] = tool_peakholdmax_clim_wrt_snr(p_ds.det_mat, [SNR_lim2, SNR_lim1], p_ds.snr_arr, p_ds.clim_arr, p_ds.N1);
  
  ## range extreme values
  v_min2 = min(p_ds.ed_mean(rr1, cc1));
  v_max2 = max(p_ds.ed_mean(rr1, cc1));
  v_min3 = min(p_ds.ed_std(rr1, cc1));
  v_max3 = max(p_ds.ed_std(rr1, cc1));
  v_min4 = min(p_ds.er_mean(rr1, cc1));
  v_max4 = max(p_ds.er_mean(rr1, cc1));
  v_min5 = min(p_ds.er_std(rr1, cc1));
  v_max5 = max(p_ds.er_std(rr1, cc1));
  
  ## plot detection state
  r_fh(1) = figure('name', 'state', 'position', [100, 100, 800, 0.62*800]);
  annotation(r_fh(1), 'textbox', [0.005, 0.005, 1, 0.1], 'fontsize', 8, ...
    'string', 'CC BY-4.0 (2025) Jakob Harden (Graz University of Technology)', 'linestyle', 'none');
  ah = axes(r_fh(1), 'tickdir', 'out');
  xlim(ah, [p_ds.clim_arr(1) - p_ds.d_clim, p_ds.clim_arr(end) + p_ds.d_clim] / p_ds.N1); 
  ylim(ah, [p_ds.snr_arr(1) - p_ds.d_snr, p_ds.snr_arr(end) + p_ds.d_snr]);
  colormap(ah, [0, 0, 0; 1, 1, 1]);
  hold(ah, 'on');
  imagesc(p_ds.clim_arr / p_ds.N1, p_ds.snr_arr, p_ds.det_mat);
  text(ah, mean([clim_lo1, clim_hi1]), mean([p_ds.snr_arr(end), SNR_lim1]), sprintf('safe choice for\nc_{lim} / N_1\n(white area)'), ...
    'horizontalalignment', 'center', 'verticalalignment', 'mid', 'fontsize', 12);
  plot(ah, [1, 1] * clim_lo1, [SNR_lim1, p_ds.snr_arr(end)], '-r', 'handlevisibility', 'off', 'linewidth', 1);
  plot(ah, [1, 1] * clim_hi1, [SNR_lim1, p_ds.snr_arr(end)], '-b', 'handlevisibility', 'off', 'linewidth', 1);
  plot(ah, [1, 1] * clim_lo2, [SNR_lim2, SNR_lim1], '-r', 'handlevisibility', 'off', 'linewidth', 1);
  plot(ah, [1, 1] * clim_hi2, [SNR_lim2, SNR_lim1], '-b', 'handlevisibility', 'off', 'linewidth', 1);
  plot(ah, [clim_lo1, clim_hi1], [1, 1] * SNR_lim1, ':k', 'handlevisibility', 'off', 'linewidth', 1);
  plot(ah, [clim_lo2, clim_hi2], [1, 1] * SNR_lim2, ':k', 'handlevisibility', 'off', 'linewidth', 1);
  plot(ah, [clim_lo1, clim_hi1], [1, 1] * p_ds.snr_arr(end), ':k', 'handlevisibility', 'off', 'linewidth', 1);
  text(ah, clim_lo1, p_ds.snr_arr(end) - 0.5, sprintf(' c_{lim} / N_1 = %.2f', clim_lo1), ...
    'fontsize', 12, 'verticalalignment', 'top');
  text(ah, clim_hi1, p_ds.snr_arr(end) - 0.5, sprintf('c_{lim} / N_1 = %.2f ', clim_hi1), ...
    'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'right');
  text(ah, mean([clim_lo1, clim_hi1]), SNR_lim1 + 1, sprintf('SNR = %.d dB', SNR_lim1), ...
    'fontsize', 12, 'verticalalignment', 'bottom', 'horizontalalignment', 'center');
  text(ah, clim_lo2, SNR_lim1 - 0.5, sprintf(' c_{lim} / N_1 = %.2f', clim_lo2), ...
    'fontsize', 12, 'verticalalignment', 'top');
  text(ah, clim_hi2, SNR_lim1 -0.5, sprintf('c_{lim} / N_1 = %.2f ', clim_hi2), ...
    'fontsize', 12, 'verticalalignment', 'top', 'horizontalalignment', 'right');
  text(ah, mean([clim_lo2, clim_hi2]), SNR_lim2 + 1, sprintf('SNR = %.d dB', SNR_lim2), ...
    'fontsize', 12, 'verticalalignment', 'bottom', 'horizontalalignment', 'center');
  hold(ah, 'off');
  colorbar(ah, 'ytick', [0, 1]);
  title(ah, sprintf('Detection state (accepted = 1 / rejected = 0)'));
  xlabel(ah, 'c_{lim} / N_1 [1]');
  ylabel(ah, 'SNR [dB]');
  
  ## differential error, mean
  r_fh(2) = figure('name', 'median', 'position', [100, 100, 800, 0.62*800]);
  annotation(r_fh(2), 'textbox', [0.005, 0.005, 1, 0.1], 'fontsize', 8, ...
    'string', 'CC BY-4.0 (2025) Jakob Harden (Graz University of Technology)', 'linestyle', 'none');
  ah = axes(r_fh(2), 'tickdir', 'out');
  xlim(ah, [p_ds.clim_arr(1) - p_ds.d_clim, p_ds.clim_arr(end) + p_ds.d_clim] / p_ds.N1);
  ylim(ah, [p_ds.snr_arr(1) - p_ds.d_snr, p_ds.snr_arr(end) + p_ds.d_snr]);
  colormap(ah, cmap2);
  hold(ah, 'on');
  imagesc(p_ds.clim_arr / p_ds.N1, p_ds.snr_arr, p_ds.ed_mean);
  text(ah, mean([clim_lo1, clim_hi1]), mean([p_ds.snr_arr(end), SNR_lim1]), sprintf('zone limits [samples]\nmin = %.2f\nmax = %.2f', v_min2, v_max2), ...
    'horizontalalignment', 'center', 'verticalalignment', 'mid', 'fontsize', 12);
  plot(ah, [1, 1] * clim_lo1, [SNR_lim1, p_ds.snr_arr(end)], ':k', 'handlevisibility', 'off', 'linewidth', 1);
  plot(ah, [1, 1] * clim_hi1, [SNR_lim1, p_ds.snr_arr(end)], ':k', 'handlevisibility', 'off', 'linewidth', 1);
  plot(ah, [clim_lo1, clim_hi1], [1, 1] * SNR_lim1, ':k', 'handlevisibility', 'off', 'linewidth', 1);
  plot(ah, [clim_lo1, clim_hi1], [1, 1] * p_ds.snr_arr(end), ':k', 'handlevisibility', 'off', 'linewidth', 1);
  hold(ah, 'off');
  colorbar();
  title(ah, sprintf('Detection error err(n_{max}) [n], mean'));
  xlabel(ah, 'c_{lim} / N_1 [1]');
  ylabel(ah, 'SNR [dB]');
  
  ## differential error, empirical standard deviation
  r_fh(3) = figure('name', 'empstd', 'position', [100, 100, 800, 0.62*800]);
  annotation(r_fh(3), 'textbox', [0.005, 0.005, 1, 0.1], 'fontsize', 8, ...
    'string', 'CC BY-4.0 (2025) Jakob Harden (Graz University of Technology)', 'linestyle', 'none');
  ah = axes(r_fh(3), 'tickdir', 'out');
  xlim(ah, [p_ds.clim_arr(1) - p_ds.d_clim, p_ds.clim_arr(end) + p_ds.d_clim] / p_ds.N1);
  ylim(ah, [p_ds.snr_arr(1) - p_ds.d_snr, p_ds.snr_arr(end) + p_ds.d_snr]);
  colormap(ah, cmap3);
  hold(ah, 'on');
  imagesc(p_ds.clim_arr / p_ds.N1, p_ds.snr_arr, p_ds.ed_std);
  text(ah, mean([clim_lo1, clim_hi1]), mean([p_ds.snr_arr(end), SNR_lim1]), sprintf('zone limits [samples]\nmin = %.2f\nmax = %.2f', v_min3, v_max3), ...
    'horizontalalignment', 'center', 'verticalalignment', 'mid', 'fontsize', 12);
  plot(ah, [1, 1] * clim_lo1, [SNR_lim1, p_ds.snr_arr(end)], ':k', 'handlevisibility', 'off', 'linewidth', 1);
  plot(ah, [1, 1] * clim_hi1, [SNR_lim1, p_ds.snr_arr(end)], ':k', 'handlevisibility', 'off', 'linewidth', 1);
  plot(ah, [clim_lo1, clim_hi1], [1, 1] * SNR_lim1, ':k', 'handlevisibility', 'off', 'linewidth', 1);
  plot(ah, [clim_lo1, clim_hi1], [1, 1] * p_ds.snr_arr(end), ':k', 'handlevisibility', 'off', 'linewidth', 1);
  hold(ah, 'off');
  colorbar();
  title(ah, sprintf('Detection error err(n_{max}) [n], emp. std.'));
  xlabel(ah, 'c_{lim} / N_1 [1]');
  ylabel(ah, 'SNR [dB]');
  
  ## relative differential error, mean
  r_fh(4) = figure('name', 'median', 'position', [100, 100, 800, 0.62*800]);
  annotation(r_fh(4), 'textbox', [0.005, 0.005, 1, 0.1], 'fontsize', 8, ...
    'string', 'CC BY-4.0 (2025) Jakob Harden (Graz University of Technology)', 'linestyle', 'none');
  ah = axes(r_fh(4), 'tickdir', 'out');
  xlim(ah, [p_ds.clim_arr(1) - p_ds.d_clim, p_ds.clim_arr(end) + p_ds.d_clim] / p_ds.N1);
  ylim(ah, [p_ds.snr_arr(1) - p_ds.d_snr, p_ds.snr_arr(end) + p_ds.d_snr]);
  colormap(ah, cmap4);
  hold(ah, 'on');
  imagesc(p_ds.clim_arr / p_ds.N1, p_ds.snr_arr, p_ds.er_mean);
  text(ah, mean([clim_lo1, clim_hi1]), mean([p_ds.snr_arr(end), SNR_lim1]), sprintf('zone limits [%%]\nmin = %.2f\nmax = %.2f', v_min4, v_max4), ...
    'horizontalalignment', 'center', 'verticalalignment', 'mid', 'fontsize', 12);
  plot(ah, [1, 1] * clim_lo1, [SNR_lim1, p_ds.snr_arr(end)], ':k', 'handlevisibility', 'off', 'linewidth', 1);
  plot(ah, [1, 1] * clim_hi1, [SNR_lim1, p_ds.snr_arr(end)], ':k', 'handlevisibility', 'off', 'linewidth', 1);
  plot(ah, [clim_lo1, clim_hi1], [1, 1] * SNR_lim1, ':k', 'handlevisibility', 'off', 'linewidth', 1);
  plot(ah, [clim_lo1, clim_hi1], [1, 1] * p_ds.snr_arr(end), ':k', 'handlevisibility', 'off', 'linewidth', 1);
  hold(ah, 'off');
  colorbar();
  title(ah, sprintf('Relative detection error err(n_{max}) / N_1 * 100 [%%], mean'));
  xlabel(ah, 'c_{lim} / N_1 [1]');
  ylabel(ah, 'SNR [dB]');
  
  ## relative differential error, empirical standard deviation
  r_fh(5) = figure('name', 'empstd', 'position', [100, 100, 800, 0.62*800]);
  annotation(r_fh(5), 'textbox', [0.005, 0.005, 1, 0.1], 'fontsize', 8, ...
    'string', 'CC BY-4.0 (2025) Jakob Harden (Graz University of Technology)', 'linestyle', 'none');
  ah = axes(r_fh(5), 'tickdir', 'out');
  xlim(ah, [p_ds.clim_arr(1) - p_ds.d_clim, p_ds.clim_arr(end) + p_ds.d_clim] / p_ds.N1);
  ylim(ah, [p_ds.snr_arr(1) - p_ds.d_snr, p_ds.snr_arr(end) + p_ds.d_snr]);
  colormap(ah, cmap5);
  hold(ah, 'on');
  imagesc(p_ds.clim_arr / p_ds.N1, p_ds.snr_arr, p_ds.er_std);
  text(ah, mean([clim_lo1, clim_hi1]), mean([p_ds.snr_arr(end), SNR_lim1]), sprintf('zone limits [%%]\nmin = %.2f\nmax = %.2f', v_min5, v_max5), ...
    'horizontalalignment', 'center', 'verticalalignment', 'mid', 'fontsize', 12);
  plot(ah, [1, 1] * clim_lo1, [SNR_lim1, p_ds.snr_arr(end)], ':k', 'handlevisibility', 'off', 'linewidth', 1);
  plot(ah, [1, 1] * clim_hi1, [SNR_lim1, p_ds.snr_arr(end)], ':k', 'handlevisibility', 'off', 'linewidth', 1);
  plot(ah, [clim_lo1, clim_hi1], [1, 1] * SNR_lim1, ':k', 'handlevisibility', 'off', 'linewidth', 1);
  plot(ah, [clim_lo1, clim_hi1], [1, 1] * p_ds.snr_arr(end), ':k', 'handlevisibility', 'off', 'linewidth', 1);
  hold(ah, 'off');
  colorbar();
  title(ah, sprintf('Relative detection error err(n_{max}) / N_1 * 100 [%%], emp. std.'));
  xlabel(ah, 'c_{lim} / N_1 [1]');
  ylabel(ah, 'SNR [dB]');
  
endfunction


function [r_clim_lo, r_clim_hi, r_ridx, r_cidx] = tool_peakholdmax_clim_wrt_snr(p_vm, p_snr, p_snr_arr, p_clim_arr, p_N1)
  ## Determine the lower and upper limit of c_lim w.r.t. given SNR limits
  ##
  ## p_vm       ... value matrix (detection error), [[<dbl>]]
  ## p_snr      ... lower and upper SNR limit [SNR_lower, SNR_upper], [<dbl>, <dbl>]
  ## p_snr_arr  ... signal-to-noise ratio variation array, [<dbl>]
  ## p_clim_arr ... c_lim variation array, [<dbl>]
  ## p_N1       ... number of samples of one full wave length, <uint>
  ## r_clim_lo  ... return: lower limit of c_lim, <dbl>
  ## r_clim_hi  ... return: upper limit of c_lim, <dbl>
  ## r_ridx     ... return: selected row index array (rows of the value matrix), [<uint>]
  ## r_cidx     ... return: selected column index array (columns of the value matrix), [<uint>]
  
  ## select block matrix from total results
  row1 = find(p_snr_arr >= p_snr(1), 1, 'first');
  row2 = find(p_snr_arr >= p_snr(2), 1, 'first');
  
  ## row index array
  r_ridx = row1 : row2;
  
  ## extract block matrix
  bm = p_vm(r_ridx, :);
  Nrow = size(bm, 1);
  
  ## assessment row vector
  ass = sum(bm);
  
  ## determine lower counter limit (x axis)
  clim_idx1 = find(ass >= Nrow, 1, 'first');
  r_clim_lo = p_clim_arr(clim_idx1) / p_N1;
  
  ## determine upper counter limit (x axis)
  clim_idx2 = find(ass >= Nrow, 1, 'last');
  r_clim_hi = p_clim_arr(clim_idx2) / p_N1;
  
  ## column index array
  r_cidx = clim_idx2 : clim_idx2;
  
endfunction


function [r_cm] = test_peakholdmax_intpcolmap(p_vm, p_cm)
  ## Interpolate colour map
  ##
  ## p_vm ... value matrix used to determine min/max, [[<dbl>]]
  ## p_cm ... colourmap matrix [Ncol x 3], [[<dbl>]]
  ## r_cm ... return: interpolated colour map matrix [64 x 3], [[<dbl>]]
  
  ## select colours from given colour map
  colhi = p_cm(end, :); # postive value color
  colcn = [1, 1, 1]; # white, center color
  collo = p_cm(1, :); # negative value color
  
  ## min/max of value matrix
  vmin = min(min(p_vm));
  vmax = max(max(p_vm));
  
  ## normlise min/max
  a = max(abs([vmin, vmax]));
  vmaxsc = vmax / a;
  vminsc = vmin / a;
  
  ## select high color
  if (abs(vmaxsc) <= 1e-6)
    c_hi = colcn;
  elseif (vmaxsc > 0)
    c_hi = colhi * vmaxsc + colcn * (1 - vmaxsc);
  else
    c_hi = collo * abs(vmaxsc) + colcn * (1 - abs(vmaxsc));
  endif
  
  ## select low color
  if (abs(vminsc) <= 1e-6)
    c_lo = colcn;
  elseif (vminsc > 0)
    c_lo = colhi * vminsc + colcn * (1 - vminsc);
  else
    c_lo = collo * abs(vminsc) + colcn * (1 - abs(vminsc));
  endif
  
  ## linear interpolation between high and low color
  N0 = 64; # number of colours
  if (vminsc <= 0) && (vmaxsc >= 0)
    #disp('dual');
    N1 = min([1, round(vmaxsc / (vmaxsc - vminsc) * N0)]);
    N2 = N0 - N1;
    cm1 = transpose(linspace(colcn, c_hi, N1 + 1))(2:end, :);
    cm2 = transpose(linspace(c_lo, colcn, N2));
    r_cm = [cm2; cm1];
  elseif (vminsc >= 0) && (vmaxsc >= 0)
    #disp('positive only');
    r_cm = transpose(linspace(c_lo, c_hi, N0));
  else
    #disp('negative only');
    r_cm = transpose(linspace(c_lo, c_hi, N0));
  endif
  
endfunction
