function [ap, prec, recall] = pascal_eval(cls, resultDir, testset, year, suffix)
% AUTORIGHTS
% -------------------------------------------------------
% Copyright (C) 2011-2012 Ross Girshick
%
% This file is part of the voc-releaseX code
% (http://people.cs.uchicago.edu/~rbg/latent/)
% and is available under the terms of an MIT-like license
% provided in COPYING. Please retain this notice and
% COPYING if you use this file (or a portion of it) in
% your project.
% -------------------------------------------------------

% cls:       object class name
% resultDir: the txt file (comp3_det_val_[cls].txt) recording the detection results and following VOC specifications

% modify the settings in voc_config.m accordingly
conf = voc_config('pascal.year', year, ...
                  'eval.test_set', testset);                
VOCopts  = conf.pascal.VOCopts;

resultFile = sprintf('%s%s%s_det_val_%s.txt', resultDir, filesep, 'comp3', cls);
if ~exist(resultFile, 'file')
    fprintf('can not find %s\n', resultFile);
    return;
end

vocResultDir = fileparts(VOCopts.detrespath);

suc = copyfile(resultFile, vocResultDir);
if suc==0
   fprintf('can not copy %s to %s\n', resultFile, vocResultDir); 
   return;
end

recall = [];
prec = [];
ap = 0;

do_eval = (str2num(year) <= 2007) | ~strcmp(testset, 'test');
if do_eval
  if str2num(year) == 2006
    [recall, prec, ap] = VOCpr(VOCopts, 'comp3', cls, true);
  else
    % Bug in VOCevaldet requires that tic has been called first
    tic;
    [recall, prec, ap] = VOCevaldet(VOCopts, 'comp3', cls, true);
  end

  % force plot limits
  ylim([0 1]);
  xlim([0 1]);

  print(gcf, '-djpeg', '-r0', [resultDir cls '_pr_' testset '_' suffix '.jpg']);
end

% save results
save([resultDir cls '_pr_' testset '_' suffix], 'recall', 'prec', 'ap');


