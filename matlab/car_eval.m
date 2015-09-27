function [rec,prec,ap] = car_eval(cls, detResultTxtFilename, gtMatFilename, minoverlap, minbboxsz, draw)
%  The PR curve is calculated based on Pascal VOC 2012 development kit

%% KITTI (eval. on the splitted half training set)
% minoverlap = 0.7;
% minbboxsz = 1000;

%% Parking_lot dataset
% minoverlap = 0.6;
% minbboxsz = 1000;

if ~exist(detResultTxtFilename, 'file') || ~exist(gtMatFilename, 'file')
    error('Can not find detection result file or groundtruth file');
end

if nargin < 4
    minoverlap = 0.5;
end

if nargin < 5
    minbboxsz = 1;
end

if nargin < 6
   draw = true; 
end

% load results
[ids,confidence,b1,b2,b3,b4]=textread(detResultTxtFilename,'%s %f %f %f %f %f');
% remove small detections
area = (b3 - b1 + 1) .* (b4 - b2 + 1);
I = area >= minbboxsz;
ids = ids(I);
BB=[b1(I) b2(I) b3(I) b4(I)]';
confidence = confidence(I);

% load gt
groundtruth = getGT(gtMatFilename, minbboxsz);

for i = 1:length(groundtruth)
    gtids{i} = groundtruth(i).name;
end

npos=0;
gt(length(gtids))=struct('BB',[],'diff',[],'det',[]);
for i=1:length(gtids)    
    gt(i).BB=groundtruth(i).BB';
    gt(i).diff = groundtruth(i).diff;
    gt(i).det = false(length(gt(i).diff),1);
    npos=npos+sum(~gt(i).diff);
end

% sort detections by decreasing confidence
[sc,si]=sort(-confidence);
ids=ids(si);
BB=BB(:,si);

% assign detections to ground truth objects
nd=length(confidence);
tp=zeros(nd,1);
fp=zeros(nd,1);
tic;
for d=1:nd
    % display progress
    if toc>1
        fprintf('%s: pr: compute: %d/%d\n',cls,d,nd);
        drawnow;
        tic;
    end
    
    % find ground truth image
    i=strmatch(ids{d},gtids,'exact');
    if isempty(i)
        error('unrecognized image "%s"',ids{d});
    elseif length(i)>1
        error('multiple image "%s"',ids{d});
    end

    % assign detection to ground truth object if any
    bb=BB(:,d);
    ovmax=-inf;
    for j=1:size(gt(i).BB,2)
        bbgt=gt(i).BB(:,j);
        bi=[max(bb(1),bbgt(1)) ; max(bb(2),bbgt(2)) ; min(bb(3),bbgt(3)) ; min(bb(4),bbgt(4))];
        iw=bi(3)-bi(1)+1;
        ih=bi(4)-bi(2)+1;
        if iw>0 & ih>0                
            % compute overlap as area of intersection / area of union
            ua=(bb(3)-bb(1)+1)*(bb(4)-bb(2)+1)+...
               (bbgt(3)-bbgt(1)+1)*(bbgt(4)-bbgt(2)+1)-...
               iw*ih;
            ov=iw*ih/ua;
            if ov>ovmax
                ovmax=ov;
                jmax=j;
            end
        end
    end
    % assign detection as true positive/don't care/false positive
    if ovmax>=minoverlap
        if ~gt(i).diff(jmax)
            if ~gt(i).det(jmax)
                tp(d)=1;            % true positive
		gt(i).det(jmax)=true;
            else
                fp(d)=1;            % false positive (multiple detection)
            end
        end
    else
        fp(d)=1;                    % false positive
    end
end

% compute precision/recall
fp=cumsum(fp);
tp=cumsum(tp);
rec=tp/npos;
prec=tp./(fp+tp);

ap=VOCap(rec,prec);

if draw
    % plot precision/recall
    plot(rec,prec,'-');
    grid;
    xlabel 'recall'
    ylabel 'precision'
    title(sprintf('class: %s, subset: %s, AP = %.3f',cls ,'test',ap));
end


% the following codes are copied from Pascal VOC 2012 devkit
function ap = VOCap(rec,prec)

mrec=[0 ; rec ; 1];
mpre=[0 ; prec ; 0];
for i=numel(mpre)-1:-1:1
    mpre(i)=max(mpre(i),mpre(i+1));
end
i=find(mrec(2:end)~=mrec(1:end-1))+1;
ap=sum((mrec(i)-mrec(i-1)).*mpre(i));


function gt = getGT(gtMatFileName, minbboxsz)
% Extract Groundtruth bounding boxes

tic;

load(gtMatFileName);

% remove images without cars
BB_old = BB;
BB = [];
testimgnames = [];
num_ids = 0;
for i = 1:length(testimagenames)
    if ~isempty(testimagenames{i})
      num_ids = num_ids + 1;
      testimgnames{num_ids} = testimagenames{i};
      BB{num_ids} = BB_old{i};
    end
end

gt=[];
for i=1:num_ids
    % display progress
    if toc>1
        fprintf('%s: pr: load: %d/%d\n','two car',i,length(testimgnames));
        drawnow;
        tic;
    end

    % read annotation
    BBox = BB{i};
    assert(~isempty(BBox));
		
    % Ground Truth
    if ~isempty(BBox)
		area = (BBox(:,4)-BBox(:,2)+1).*(BBox(:,3)-BBox(:,1)+1);
		% we set small bbox (lower than minbboxsz) as difficult examples.
		Diff = zeros(1, size(BBox,1));
		I = find(area < minbboxsz);
		Diff(I) = 1;
        gt(i).BB=BBox;
        gt(i).diff = Diff;
        gt(i).im = testimgnames{i};
        gt(i).name = testimgnames{i};
    end
end