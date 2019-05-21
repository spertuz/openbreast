function dataset = readBCDR(csvfile)
% Read image date from BCDR dataset
% Sintax:
%     dataset = readBCDR(csvfile)
% Inputs:
%     csvfile,    path to *.csv file with image data
% Outputs:
%     dataset,    a table with data from each image
%                 in the dataset.
%                 
%  S. Pertuz
%  Jan09/2018

% determine whether the input is a _img or _outfile file:
[~, fname] = fileparts(csvfile);
k = strfind(fname, '_');
isImg = strcmp('_img', fname(k(end):end));

if isImg
    %read img file
    f = fopen(csvfile);
    header = textscan(f, '%s', 8, 'delimiter', ',');
    lines = textscan(f, '%d %d %d %s %s %d %d %d', 'delimiter', ','); 
    fclose(f);
else
    %read outline file:
    f = fopen(csvfile);
    header = textscan(f, '%s', 19, 'delimiter', ',');
    lines = textscan(f, '%d %d %d %d %d %d %s %s %s %d %d %d %d %d %d %d %d %d %s', ...
        'delimiter', ',');
    fclose(f);
end

data1 = cell2struct(lines, header{1}, 2); 

fnames = {'patient_id','study_id','series'};
for n = 1:length(fnames)
    if isfield(data1, fnames{n})
        dataset.(fnames{n}) = data1.(fnames{n});    
    end
end
%breast density
if isfield(data1, 'density')
    dataset.density = data1.density;
else
    dataset.density = zeros(length(dataset.patient_id), 1, 'uint8');
end

%lession age
if isfield(data1, 'age')
    dataset.age = data1.age;
else
    dataset.age = nan(length(dataset.patient_id), 1);
end

%lession id
if isfield(data1, 'lession_id')
    dataset.lession_id = data1.lession_id;
else
    dataset.lession_id = cell(length(dataset.patient_id), 1);
end

%correct image path's
if isfield(data1, 'image_filename')
    [srcpath, ~, ~] = fileparts(csvfile);
    fnames = data1.image_filename;
    dataset.path = cell(length(fnames), 1);
    for n = 1:length(fnames)
        dataset.path{n} = [srcpath, filesep, fnames{n}];
    end
end
    

%correct mammo class
if isfield(data1, 'classification')
    mclass = data1.classification;
    dataset.class = false(length(mclass), 1);
    for n = 1:length(mclass)
        dataset.class(n) = strcmp(mclass{n}, 'Malign');
    end
else
    dataset.class = false(size(dataset.patient_id));
end

%correct image view:
if isfield(data1, 'image_view')
    view = data1.image_view;
    view_list = {'RCC', 'LCC', 'RMLO', 'LMLO'};
    dataset.view = cell(length(view), 1);
    for n = 1:length(view)
        dataset.view{n} = view_list{view(n)};
    end
elseif isfield(data1, 'image_type_name')
    view = data1.image_type_name;
    view_list = {'RCC', 'LCC', 'RMLO', 'LMLO'};
    view_code = {'RCC', 'LCC', 'RO', 'LO'};
    dataset.view = cell(length(view), 1);
    for n = 1:length(view)
        dataset.view{n} = view_list{strcmp(view_code, view{n})};
    end
end

%correct x-points
if isfield(data1, 'lw_x_points')
    xpts = data1.lw_x_points;
    ypts = data1.lw_y_points;
    
    dataset.xpoints = cell(length(xpts), 1);
    dataset.ypoints = cell(length(ypts), 1);
    for n = 1:length(xpts)
        dataset.xpoints{n} = str2num(xpts{n});
        dataset.ypoints{n} = str2num(ypts{n});
    end
else
    dataset.xpoints = cell(length(dataset.patient_id), 1);
    dataset.ypoints = cell(length(dataset.patient_id), 1);
end

dataset = struct2table(dataset);