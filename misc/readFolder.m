function imdata = readFolder(srcdir, tflag, src)
% Read metadata from Measurement Challenge directory
% Sintax:
%     imdata = readFolder(srcdir)
%     imdata = readFolder(srcdir, tflag)
% Inputs:
%     srcdir,     path to source directory
%     tflag,      boolean flag. True (default) for training set.
% Outputs:
%     imdata,     a table with information of each image in the
%                 dataset. The structure has the following fields:
%         ID,             patient ID
%         year1,          year of birth
%         year2,          year of mammo
%         ethnicity,      1 (Asian subcontinent), 2 (European), 3 (hispanic),
%                         4 (missing or other), 5 (South-east asian), 6
%                         (african).
%         bmi,            Body mass index
%         side,           'L' or 'R'
%         view,           'MLO' or 'CC'
%         isRaw,          boolean flag. 1 for RAW, 0 for PROCESSED
%         vendor,         1 (Afga), 2 (fujifilm), 3 (hologic), 4 (konica), 5 (philips), 
%                         6 (siemens), 7 (sectra) or 8 (GE). 
%         impath,         path to dicom image
%         class,          0 for control, 1 for 'invasive', 2 for 'in situ'
%         isTrain,        boolean flag. true for training images.
%         source,         1 (australia), 2 (malasya), 3 (norway), 4 (UK), 
%                         5 (USA).
% 
% 
% S. Pertuz
% Sep27/2017

if nargin<2
    tflag = true;
end

if tflag
    t = readtable([srcdir, filesep, 'CoreDataTraining.csv']);
    imdir = [srcdir, filesep, 'Images', filesep, 'Training', filesep];
    isTrain = table(true(size(t,1), 1), 'VariableNames', {'isTrain'});
    source = table(src*ones(size(t,1), 1), 'VariableNames', {'source'});  
    imdata = t(:, [1 2 3 7 10 11 12 13 14 15 8]);
    imdata = [imdata, isTrain, source];
    imdata.Properties.VariableNames = {'ID', 'year1', 'year2', 'ethnicity',...
       'BMI', 'side', 'view', 'isRaw', 'vendor', 'path',  'class', 'isTrain', 'source'};
else
    t = readtable([srcdir, filesep, 'CoreDataTesting.csv']);
    imdir = [srcdir, filesep, 'Images', filesep, 'Testing', filesep];
    isTrain = table(false(size(t,1), 1), 'VariableNames', {'isTrain'});
    source = table(src*ones(size(t,1), 1), 'VariableNames', {'source'});
    class = table(nan(size(t,1), 1), 'VariableNames', {'class'});
    imdata = [t(:, [1 2 3 6 7 8 9 10 11 12]), class, isTrain, source];
    imdata.Properties.VariableNames = {'ID', 'year1', 'year2', 'ethnicity',...
       'BMI', 'side', 'view', 'isRaw', 'vendor', 'path',  'class', 'isTrain', 'source'};
end


% Retrieve list of ethnicities, manufacturers
% and diagnosis:
%Ethnicity list
elist = unique(imdata.ethnicity(:));
e_class = zeros(length(elist), 1);
fprintf('Ethnicities:\n')
for n = 1:length(elist)
    e_class(n) = input(sprintf('%s: ', elist{n}));
end

%Manufacturer list
mlist = unique(imdata.vendor(:));
m_class = zeros(length(mlist), 1);
fprintf('\nManufacturers:\n')
for n = 1:length(mlist)
    m_class(n) = input(sprintf('%s: ', mlist{n}));
end

% Diagnosis list
if tflag
    slist = unique(imdata.class(:));
    s_class = zeros(length(slist), 1);
    fprintf('\nDianosis:\n')
    for n = 1:length(slist)
        s_class(n) = input(sprintf('%s: ', slist{n}));
    end
end

for n = 1:size(imdata, 1)
    %Fix mammo side:
    if strcmpi(imdata.side{n}, 'Left')
        imdata.side{n} = 'L';
    else
        imdata.side{n} = 'R';
    end
    
    %Fix mammo type:
    imdata.isRaw{n} = ~strcmpi(imdata.isRaw{n}, 'Processed');
    
    %Fix file ID:
    fid = imdata.path{n};
    i = strfind(fid, '/');
    if ~isempty(i)
        fid = fid(i(end)+1:end);
    end
    
    %Fix path
    [folder, ~] = strtok(fid, '_');
    impath = [imdir, folder, filesep, fid];
    [fpath, fname, ~] = fileparts(impath);
    imdata.path{n} = [fpath, filesep, fname, '.dcm'];
    
    %Diagnosis
    if tflag
        imdata.class{n} = s_class(strcmp(imdata.class{n}, slist));
    end
    
    %Ethnicity
    imdata.ethnicity{n} = e_class(strcmp(imdata.ethnicity{n}, elist));
    
    %Manufacturer
    imdata.vendor{n} = m_class(strcmp(imdata.vendor{n}, mlist));
end

% Convert output to structure:
imdata = table2struct(imdata, 'ToScalar', true);
imdata.ethnicity = cell2mat(imdata.ethnicity);
imdata.isRaw = cell2mat(imdata.isRaw);
imdata.vendor = cell2mat(imdata.vendor);
if tflag
    imdata.class = cell2mat(imdata.class);
end

%Display summary:
%RAW images:
select = imdata.isRaw(:);
fprintf('\nRAW images: %d ******\n', sum(select))
age = imdata.year2(select)-imdata.year1(select);
fprintf('Age:       [%d-%d] (%2.0d)\n', min(age), max(age), median(age))
eth = imdata.ethnicity(select);
fprintf('Ethnicity: 1:%d, 2:%d, 3:%d, 4:%d, 5:%d, 6:%d\n',...
    sum(eth==1), sum(eth==2), sum(eth==3), sum(eth==4), sum(eth==5), sum(eth==6));
bmi = imdata.BMI(select);
fprintf('BMI:       [%2.2f-%2.2f] (%2.2f)\n', min(bmi), max(bmi), mean(bmi))
v = imdata.vendor(select);
fprintf('Vendor:    1:%d, 2:%d, 3:%d, 4:%d, 5:%d, 6:%d, 7:%d, 8:%d\n',...
    sum(v==1), sum(v==2), sum(v==3), sum(v==4), sum(v==5), sum(v==6), sum(v==7), sum(v==8));
if tflag
    d = imdata.class(select);
    fprintf('Diagnosis: control: %d, invasive: %d, in situ: %d\n', sum(d==0), sum(d==1), sum(d==2));
end

%PROCESSED images:
select = ~imdata.isRaw(:);
fprintf('\nPROCESSED images: %d ******\n', sum(select))
age = imdata.year2(select)-imdata.year1(select);
fprintf('Age:       [%d-%d] (%2.0d)\n', min(age), max(age), median(age))
eth = imdata.ethnicity(select);
fprintf('Ethnicity: 1:%d, 2:%d, 3:%d, 4:%d, 5:%d, 6:%d\n',...
    sum(eth==1), sum(eth==2), sum(eth==3), sum(eth==4), sum(eth==5), sum(eth==6));
bmi = imdata.BMI(select);
fprintf('BMI:       [%2.2f-%2.2f] (%2.2f)\n', min(bmi), max(bmi), mean(bmi))
v = imdata.vendor(select);
fprintf('Vendor:    1:%d, 2:%d, 3:%d, 4:%d, 5:%d, 6:%d, 7:%d, 8:%d\n',...
    sum(v==1), sum(v==2), sum(v==3), sum(v==4), sum(v==5), sum(v==6), sum(v==7), sum(v==8));
if tflag
    d = imdata.class(select);
    fprintf('Diagnosis: control: %d, invasive: %d, in situ: %d\n', sum(d==0), sum(d==1), sum(d==2));
end

imdata = struct2table(imdata);
end