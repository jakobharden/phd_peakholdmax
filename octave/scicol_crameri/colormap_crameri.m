## Load colormaps desgined by Fabio Crameri et. al., version 8
##
## Usage: [r_cm, r_fh] = colormap_crameri(p_cmn, p_num)
##
## p_cmn ... colourmap name, <str>
## p_num ... number of colours, <uint>
## r_cm  ... return: colourmap [Ncol, 3], [[<dbl>]]
## r_fh  ... return: figure handle, colourmap preview, optional, <uint>
##
## Note: If p_cmn = 'list', print available colour map names on screen.
##
## Website: https://www.fabiocrameri.ch/colourmaps/
## Source: https://zenodo.org/records/8409685
## All colormap files by Fabio Crameri et. al. available at Zenodo are licensed under the MIT licence.
##
#######################################################################################################################
## LICENSE
##
##    Copyright (C) 2024 Jakob Harden (jakob.harden@tugraz.at, Graz University of Technology, Graz, Austria)
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
## This file is part of the PhD thesis of Jakob Harden.
##
function [r_cm, r_fh] = colormap_crameri(p_cmn, p_num)
  
  ## check arguments
  if (nargin < 2)
    p_num = [];
  endif
  if (nargin < 1)
    p_cmn = [];
  endif
  if isempty(p_num)
    p_num = 16;
  endif
  if isempty(p_cmn)
    p_cmn = 'batlow';
  endif
  p_num = min([max([p_num, 1]), 256]); # limit number
  
  ## print available color map names
  switch (p_cmn)
    case 'list'
      ## list available color maps
      disp('Available color maps (Scientific Colors by Fabio Crameri et. al., Version 8)');
      disp('  acton');
      disp('  bamako');
      disp('  bam');
      disp('  bamO');
      disp('  batlowK');
      disp('  batlow');
      disp('  batlowW');
      disp('  berlin');
      disp('  bilbao');
      disp('  broc');
      disp('  brocO');
      disp('  buda');
      disp('  bukavu');
      disp('  cork');
      disp('  corkO');
      disp('  davos');
      disp('  devon');
      disp('  fes');
      disp('  glasgow');
      disp('  grayC');
      disp('  hawaii');
      disp('  imola');
      disp('  lajolla');
      disp('  lapaz');
      disp('  lipari');
      disp('  lisbon');
      disp('  managua');
      disp('  navia');
      disp('  naviaW');
      disp('  nuuk');
      disp('  oleron');
      disp('  oslo');
      disp('  roma');
      disp('  romaO');
      disp('  tofino');
      disp('  tokyo');
      disp('  turku');
      disp('  vanimo');
      disp('  vikO');
      return;
  endswitch
  
  ## file name, file path
  colfn = sprintf('%s.mat', p_cmn);
  colfp = fullfile('.', 'scicol_crameri', colfn);
  
  ## check if file exists
  if not(exist(colfp, 'file') == 2)
    error('Color map "%s" undefined!', p_cmn);
  endif
  
  ## load colormap
  dat = load(colfp, '-v6');
  cm = getfield(dat, p_cmn);
  
  ## interpolate colors
  xi = linspace(1, size(cm, 1), p_num)';
  xs = [1 : 256]';
  r_cm = zeros(p_num, 3);
  r_cm(:, 1) = interp1(xs, cm(:, 1), xi, 'linear');
  r_cm(:, 2) = interp1(xs, cm(:, 2), xi, 'linear');
  r_cm(:, 3) = interp1(xs, cm(:, 3), xi, 'linear');
  
  ## plot preview
  if isargout(2)
    r_fh = figure('name', 'colourmap_crameri', 'position', [300, 300, 640, 640 / 1.62]);
    ah = axes(r_fh, 'position', [0.05, 0.1, 0.9, 0.8], 'clipping', 'off');
    axis(ah, 'off');
    xlim(ah, [0.5, p_num + 0.5]);
    hold(ah, 'on');
    for j = 1 : p_num
      fill(ah, [j - 0.5, j + 0.5, j + 0.5, j - 0.5], [0, 0, 1, 1], r_cm(j, :), 'edgecolor', 'none');
    endfor
    rectangle(ah, 'Position', [0.5, 0, p_num, 1], 'edgecolor', 'k', 'facecolor', 'none');
    hold(ah, 'off');
    title(ah, sprintf('Scientific colour map: %s (by Fabio Crameri et. al.)', p_cmn), 'fontsize', 12);
    annotation('textbox', [0.05, 0.01, 0.9, 0.05], 'string', 'CC BY-4.0 Jakob Harden (Graz University of Technology), 2024', 'linestyle', 'none');
  endif
  
endfunction
