function [GLRLMS, SI] = grayrlmatrix(I, Offset, NL, GL, mask)
% from: https://se.mathworks.com/matlabcentral/fileexchange/17482-gray-level-run-length-matrix-toolbox
%  Description
%  -------------------------------------------
%   Computes the graylevel run length (GLRL) matrix used for texture
%   analysis of an image using zigzag scan method.
%
%   [GLRLMS,SI]= grayrlmatrix(I, Offset, NL, GL, mask) returns one or more
%   gray-level run-length matrices, depending on the values of the optional
%   parameter/value pairs.
%
%  ------------------------------------------
%  Example
%  ------------------------------------------
% I =[1     1     1     2     2
%      3     4     2     2     3
%      4     4     4     4     4
%      5     5     3     3     3
%      1     1     3     4     5]
% [GLRLMS,SI] = grayrlmatrix(I,'NumLevels',5,'G',[])
% I =
%      1     1     1     2     2
%      3     4     2     2     3
%      4     4     4     4     4
%      5     5     3     3     3
%      1     1     3     4     5
% GLRLMS(:,:,1) =
%      0     1     1     0     0
%      0     2     0     0     0
%      3     0     1     0     0
%      2     0     0     0     1
%      1     1     0     0     0
% GLRLMS(:,:,2) =
%      5     0     0     0     0
%      0     2     0     0     0
%      4     1     0     0     0
%      5     1     0     0     0
%      3     0     0     0     0
% GLRLMS(:,:,3) =
%      5     0     0     0     0
%      2     1     0     0     0
%      4     1     0     0     0
%      5     1     0     0     0
%      3     0     0     0     0
% GLRLMS(:,:,4) =
%      5     0     0     0     0
%      4     0     0     0     0
%      6     0     0     0     0
%      5     1     0     0     0
%      3     0     0     0     0
% SI =
%      1     1     1     2     2
%      3     4     2     2     3
%      4     4     4     4     4
%      5     5     3     3     3
%      1     1     3     4     5
% -------------------------------------------
% See also zigzag rle_0 rle_45
% -------------------------------------------
% Author:
% -------------------------------------------
%    (C)Xunkai Wei <xunkai.wei@gmail.com>
%    Beijing Aeronautical Technology Research Center
%    Beijing %9203-12,100076
% -------------------------------------------
% History:
% -------------------------------------------
% Creation: beta  Date: 01/10/2007
% Revision: 1.0   Date: 14/11/2007
% Revision: 2.0   Date: 20/01/2016
% -------------------------------------------
% Bug Fixed:
% -------------------------------------------
% 1.Issue wrong results for nonsquare matrix,now output cells instead of
%   multi-dim arrays
% 2.Add support for inputs checking inspired by MATLAB style
% 3.It includes mask binary input specifying pixels of interest 
% 

% Scale I so that it contains integers between 1 and NL.
if GL(2) == GL(1)
    SI = ones(size(I));
else
    slope = (NL-1) / (GL(2) - GL(1));
    intercept = 1 - (slope*(GL(1)));
    SI = round(imlincomb(slope,I,intercept,'double'));
end

% Clip values if user had a value that is outside of the range, e.g., double
% image = [0 .5 2;0 1 1]; 2 is outside of [0,1]. The order of the following
% lines matters in the event that NL = 0.
SI(SI > NL) = NL;
SI(SI < 1) = 1;

% total numbers of directions
numOffsets = length(Offset);
GLRLMS = cell(numOffsets, 1);

if NL ~= 0
    % make direction matrix for all given directions
    for k = 1 : numOffsets
        GLRLMS{k} = computeGLRLM(SI,Offset(k),NL,mask);
    end
else
    GLRLMS = [];
end
end

% --------------------------------------------------------------------
function oneGLRLM = computeGLRLM(si,offset,nl,mask)
% For given direction, compute the run length matrix
switch offset
    case 1
        % 0 degree
        oneGLRLM = rle_0(si,nl,mask);
    case 2
        % 45 degree
        seq = zigzag(si);
        seq_mask = zigzag(mask);
        oneGLRLM  = rle_45(seq,nl,seq_mask);
    case 3
        % 90 degree
        oneGLRLM = rle_0(si',nl,mask');
    case 4
        % 135 degree
        seq = zigzag(fliplr(si));
        seq_mask = zigzag(fliplr(mask));
        oneGLRLM = rle_45(seq,nl,seq_mask);
    otherwise
        error('Only 4 directions supported')
end
end

% ---------------------------------------------
function oneglrlm = rle_0(si,NL,mask)
% RLE   image gray level Run Length matrix for 0degree
%    
% Author:
% ---------------------------------------------
%    (C)Xunkai Wei <xunkai.wei@gmail.com>
%    Beijing Aeronautical Technology Research Center
%    Beijing %9203-12,10076
%
%    Germ�n F. Torres V.    (Revision 2.0)
%    Universidad Industrial de Santander
%    Bucaramanga
%
% History:
% ---------------------------------------------
% Creation: beta  Date: 01/11/2007 
% Revision: 1.0   Date: 10/11/2007
% Revision: 2.0   Date: 20/01/2016
%   It includes mask binary input specifying pixels of interest 
% ---------------------------------------------


% Assure row number is exactly the gray level
[m,n]=size(si);

oneglrlm=zeros(NL,n);

for i=1:m
    x=si(i,:);
    x = x(mask(i,:));
    if ~isempty(x)
        % run length Encode of each vector
        index = [ find(x(1:end-1) ~= x(2:end)), length(x) ];
        len = diff([ 0 index ]); % run lengths
        val = x(index);          % run values
        temp =accumarray([val;len]',1,[NL n]);% compute current numbers (or contribution) for each bin in GLRLM
        oneglrlm = temp + oneglrlm; % accumulate each contribution
    end
end
end

function oneglrlm = rle_45(seq,NL,seq_mask)
% RLE   image gray level Run Length matrix for 45 and 135
% This file is to handle the zigzag scanned sequence for 45 or 135 degree
% direction. Note for 135, just swap the left and the right colum
% Author:
% ---------------------------------------------
%    (C)Xunkai Wei <xunkai.wei@gmail.com>
%    Beijing Aeronautical Technology Research Center
%    Beijing %9203-12,10076
%
%    Germ�n F. Torres V.    (Revision 2.0)
%    Universidad Industrial de Santander
%    Bucaramanga
%
% History:
% ---------------------------------------------
% Creation: beta  Date: 01/11/2007 
% Revision: 1.0   Date: 10/11/2007
% Revision: 2.0   Date: 20/01/2016
%   It includes sequence of mask binary input specifying pixels of interest 
% ---------------------------------------------


% Assure row number is exactly the gray leve;
% number of seqence
m =length(seq);
% number to store the possible max coloums
n = findmaxnum(seq);
%

oneglrlm=zeros(NL,n);

for i=1:m
    x = seq{i};
    if isscalar(x)
        if seq_mask{i}
                index = [ find(x(1:end-1) ~= x(2:end)), length(x) ];
                len = diff([ 0 index ]); % run lengths
                val = x(index);          % run values
                temp =accumarray([val;len]',1,[NL n]);% compute current numbers (or contribution) for each bin in GLRLM
                oneglrlm = temp + oneglrlm; % accumulate each contribution
        end
    else
        x=x(logical(seq_mask{i}));
        if ~isempty(x)
            % run length Encode of each vector
            index = [ find(x(1:end-1) ~= x(2:end)), length(x) ];
            len = diff([ 0 index ]); % run lengths
            val = x(index);          % run values
            temp =accumarray([val;len]',1,[NL n]);% compute current numbers (or contribution) for each bin in GLRLM
            oneglrlm = temp + oneglrlm; % accumulate each contribution
        end
    end
end
end

% ---------------------------------------------
function seq = zigzag(SI)
%
%  Description:
%  ------------
%  This function is used to build the corresponding sequences of a given
%  scaled gray level image matrix from 45' degree direction. The whole process is using zigzag method
%  It can handle nonsquare image matrix
%
% Author:
% -------
%    (C)Xunkai Wei <xunkai.wei@gmail.com>
%    Beijing Aeronautical Technology Research Center
%    Beijing %9203-12,10076
%
% History:
%  -------
% Creation: beta  Date: 01/11/2007
% Revision: 1.0   Date: 12/11/2007
% 
% Trick: all the sequence starts or ends lie on the boundary.

% initializing the variables
%----------------------------------
c = 1; % initialize colum indicator
r = 1; % initialize row   indicator

rmin = 1; % row   boundary checker
cmin = 1; % colum boundary checker

rmax = size(SI, 1); % get row numbers
cmax = size(SI, 2); % get colum numbers

%
i = 1; % counter for current ith element
j = 1; % indicator for determining sequence interval

% intialize sequence mark
sq_up_begin=1;

sq_down_begin=1;

% % Output contain value and its flag status
% the first row contain value
% the second row contain its flag
output = zeros(1, rmax * cmax);
% sequence counter
%
% % Position Matrix
% position =zeros(1, rmax * cmax);
%----------------------------------

while ((r <= rmax) && (c <= cmax))

    % for current point, judge its zigzag direction up 45, or down 45, or
    % 0,or down 90

    if (mod(c + r, 2) == 0)      % up 45 direction
        %  if we currently walk to the left first colum
        if (r == rmin)
            % First, record current point
            output(i) = SI(r, c);
            % if we walk to right last colum
            if (c == cmax)
                % add row number move straight down 90
                r   = r + 1;
                sq_up_end = i;
                sq_down_begin = i+1;
                seq{j}=output(sq_up_begin:sq_up_end);
                j = j + 1;
                %

            else
                % Continue to move to next (1,c+1) point
                % This next point should be the begin point of next sequence
                c = c + 1;
                sq_up_end = i;
                sq_down_begin = i+1;

                seq{j}=output(sq_up_begin:sq_up_end);

                j = j + 1;

            end;

            % add couter
            i = i + 1;
            % if we currently walk to the last column
        elseif ((c == cmax) && (r < rmax))
            % first record the point
            output(i) = SI(r, c);
            % then move straight down to next row
            r = r + 1;
            
            sq_up_end = i;
            seq{j}=output(sq_up_begin:sq_up_end);
            sq_down_begin =i+1;
            j=j+1;
                        
            % add counter
            i = i + 1;
            % all other cases i.e. nonboundary points
        elseif ((r > rmin) && (c < cmax))
            output(i) = SI(r, c);
            % move to next up 45 point
            r = r - 1;
            c = c + 1;
            % add counter
            i = i + 1;
        end;
        % down 45 direction
    else
        % if we walk to the last row
        if ((r == rmax) && (c <= cmax))
            % firstly record current point
            output(i) = SI(r, c);
            % move right to next point
            c = c + 1;
            sq_down_end = i;
            seq{j}=output(sq_down_begin:sq_down_end);
            sq_up_begin =i+1;
            j = j + 1;
            % add counter
            i = i + 1;
            % if we walk to the first column
        elseif (c == cmin)
            %first record current point
            output(i) = SI(r, c);
            %
            if (r == rmax)
                c = c + 1;
                
                sq_down_end = i;
                seq{j}=output(sq_down_begin:sq_down_end);
                sq_up_begin =i+1;
                j = j + 1;

            else
                r = r + 1;
                % record sequence end
                sq_down_end = i;
                seq{j}=output(sq_down_begin:sq_down_end);
                sq_up_begin =i+1;
                j = j + 1;

            end;

            i = i + 1;
            % all other cases without boundary point
        elseif ((r < rmax) && (c > cmin))
            %
            output(i) = SI(r, c);
            %             position(i) = sub2ind(SI,r,c);
            r = r + 1;
            c = c - 1;
            % keep down_info
            i = i + 1;
        end;

    end;

    if ((r == rmax) && (c == cmax))          % bottom right element
        output(i) = SI(r, c);
        sq_end = i;
        seq{j}=output(sq_end);
        %         position(i) = sub2ind(SI,r,c);
        break
    end;


end;
end

function maxnum=findmaxnum(seq)
%  This function obtain the maximum numbers of the given sequence
%  note the sequence is stored in cell mode
%
%
% See also zigzag
%
% Author:
% ---------------------------------------------
%    (C)Xunkai Wei <xunkai.wei@gmail.com>
%    Beijing Aeronautical Technology Research Center
%    Beijing %9203-12,10076
% History:
% ---------------------------------------------
% Creation: beta  Date: 01/10/2007
% Revision: 1.0   Date: 12/11/2007
%
if iscell(seq)

    numseq=length(seq);
    maxnum=1;
    for i=1:numseq
        temp = seq{i};
        numseq = length(temp);
        if numseq > maxnum
            maxnum =numseq;
        end
    end
else
    error('I was only designed to handle cell sequence')
end
end