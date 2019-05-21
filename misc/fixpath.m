function fpath = fixpath(srcpath, prefix, no, ext)
% Fix image path
% Sintax:
%     fpath = fixpath(srcpath, prefix, no)
%     fpath = fixpath(srcpath, prefix, no, ext)
% Inputs:
%     srcpath,    input path
%     prefix,     prefix for replacing in path
%     no,         index from which the srcpath is kept
%     ext,        file extension (default dcm).
% 
% S. Pertuz
% Nov17/2017

if nargin<4
    ext = '.dcm';
end

if iscell(srcpath)
    fpath = cell(size(srcpath));
    for n = 1:length(srcpath)
        fpath{n} = fixpath(srcpath{n}, prefix, no, ext);
    end
    return
end

fpath = [prefix, filesep, srcpath(no:end)];
[srcdir, fname, ~] = fileparts(fpath);
fpath = [srcdir, filesep, fname, ext];