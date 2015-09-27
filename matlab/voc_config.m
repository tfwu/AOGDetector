function conf = voc_config(varargin)
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

%
% ~~~~~~~~~~~~~~~~~~~~~~ BASIC SETUP ~~~~~~~~~~~~~~~~~~~~~~
% Please read the next few lines

% Parent directory that everything (model cache, VOCdevkit) is under
% I recommend making a symlink to your BASE_DIR named 'cachedir'.
% e.g., cachedir -> /var/tmp/rbg/
BASE_DIR    = '/home/tfwu/Data/'; 

% PASCAL dataset year to use
PASCAL_YEAR = '2007';

% The code will look for your PASCAL VOC devkit in
% BASE_DIR/VOC<PASCAL_YEAR>/VOCdevkit
% e.g., /var/tmp/rbg/VOC2007/VOCdevkit
% If you have the devkit installed elsewhere, you may want to
% create a symbolic link.

% Configuration structure
conf = [];

% Clobber with overrides passed in as arguments
for i = 1:2:length(varargin)
    key = varargin{i};
    val = varargin{i+1};
    eval(['conf.' key ' = val;']);
end

% -------------------------------------------------------------------
% PASCAL VOC configuration
% -------------------------------------------------------------------
% Parent directory that everything (model cache, VOCdevkit) is under
conf = cv(conf, 'paths.base_dir', BASE_DIR);

% Configure the PASCAL VOC dataset year
conf = cv(conf, 'pascal.year', PASCAL_YEAR);

% Directory with PASCAL VOC development kit and dataset
conf = cv(conf, 'pascal.dev_kit', [conf.paths.base_dir '/VOC' ...
    conf.pascal.year '/VOCdevkit/']);

% For INRIA person
% conf = cv(conf, 'pascal.dev_kit', [conf.paths.base_dir '/INRIA_PASCAL/VOCdevkit/']);

if exist(conf.pascal.dev_kit, 'dir') == 0   
    msg = sprintf(['~~~~~~~~~~~ Hello ~~~~~~~~~~~\n' ...
        'voc-release5 is not yet configured for learning. \n' ...
        'You can still run demo.m, but please read \n' ...
        'the section "Using the learning code" in README. \n' ...
        '(Could not find the PASCAL VOC devkit in %s)'], ...
        conf.pascal.dev_kit);
    fprintf([msg '\n\n']);   
    return;
end

% VOCinit brings VOCopts into scope
conf.pascal.VOCopts = get_voc_opts(conf);
conf.pascal.VOCopts.testset = conf.eval.test_set;

% -------------------------------------------------------------------
% Returns the 'VOCopts' variable from the VOCdevkit. The path to the
% devkit is also added to the matlab path.
function VOCopts = get_voc_opts(conf)
% cache VOCopts from VOCinit
persistent voc_opts;

key = conf.pascal.year;
if isempty(voc_opts) || ~voc_opts.isKey(key)
    if isempty(voc_opts)
        voc_opts = containers.Map();
    end
    tmp = pwd;
    cd(conf.pascal.dev_kit);
    addpath([cd '/VOCcode']);
    VOCinit;
    cd(tmp);
    voc_opts(key) = VOCopts;
end
VOCopts = voc_opts(key);


% -------------------------------------------------------------------
% Does nothing if conf.key exists, otherwise sets conf.key to val
function conf = cv(conf, key, val)
try
    eval(['conf.' key ';']);
catch
    eval(['conf.' key ' = val;']);
end

