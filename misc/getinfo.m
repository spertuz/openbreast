function info = getinfo(dicom_image)
% Retrieve information DICOM file
% SINTAX:
%     INFO = getinfo(DICOM_IMAGE);
% Retrieves and prints the information about the dicom image found
% in the path DICOM_IMAGE. The output is a structure with the 
% following fields:
%     info.ID,        An integer with the accession number
%     info.istomo,    A flag. True for tomo images.
%     info.iscc,      A flag. True for CC views.
%     info.view,      A string. {CC, MLO}
%     info.side,      A string. {L or R}
%     info.israw,     A flag. True for raw images.
%     info.tomo,      A string. 'RC' (reconstruction)  or 'PR' (projection)
%     info.source,    A string with the manufacturer
%     info.target,    A string with the target material 'MO' or 'RH'
%     info.filter,    A string with the filter material 'MO' or 'RH'
%     info.KVP,       An integer with the KVP
%     info.psize,     A float with the pixel size (mm)
%     info.H,         An integer with the breast thickness (in mm)
%     info.cforce,    Compression force
%     info.exposure,  Exposure time
%     info.dose,      Organ dose
%  
% Said Pertuz
% Apr18/2014

info.ID = [];
info.source = [];
info.israw = false;
info.view = [];
info.side = [];
info.tomo = [];
info.istomo = false;
info.ismlo = false;
info.iscc = false;

if ~exist(dicom_image, 'file')
    warning('the file does not exist!')
    return
end

full_info = dicominfo(dicom_image, 'UseDictionaryVR', true);

%ISTOMO
if ~isfield(full_info, 'SeriesDescription')||isempty(strfind(full_info.SeriesDescription, 'Tomosynthesis'))
    info.istomo = false;
    info.tomo = [];
else
    info.istomo = true;
    if strfind(full_info.SeriesDescription, 'Reconstruction')
        info.tomo = 'RC';
    else
        info.tomo = 'PR';
    end
    info.israw = strfind(full_info.SeriesDescription, 'Raw');        
end

%ISRAW
if isfield(full_info, 'PresentationIntentType')
    if strcmp(full_info.PresentationIntentType,'FOR PRESENTATION')
        info.israw = false;
    else
        info.israw = true;
    end
end

%ACCESION NUMBER
if isfield(full_info, 'AccessionNumber')
    info.ID = full_info.AccessionNumber;
end

%VIEW: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(full_info, 'ViewPosition')&&strcmp(full_info.ViewPosition, 'MLO')
    info.view = full_info.ViewPosition;    
elseif isfield(full_info, 'ViewPosition')&&strcmp(full_info.ViewPosition, 'CC')
    info.view = 'CC';    
elseif isfield(full_info, 'SeriesDescription')&&contains(full_info.SeriesDescription, 'CC')
    info.view = 'CC';    
elseif isfield(full_info, 'SeriesDescription')&&contains(full_info.SeriesDescription, 'MLO')
    info.view = 'MLO';    
elseif isfield(full_info, 'ProtocolName')&&~isempty(full_info.ProtocolName)&&contains(full_info.ProtocolName, 'MLO')
    info.view = 'MLO';    
elseif isfield(full_info, 'ProtocolName')&&~isempty(full_info.ProtocolName)&&contains(full_info.ProtocolName, 'CC')
    info.view = 'CC';    
elseif isfield(full_info, 'ViewPosition')&&~isempty(full_info.ViewPosition)
    info.view = full_info.ViewPosition;    
elseif isfield(full_info, 'AcquisitionDeviceProcessingDescription')&&~isempty(full_info.AcquisitionDeviceProcessingDescription)
    tmp = full_info.AcquisitionDeviceProcessingDescription;
    if contains(tmp, 'CC')
        info.view = 'CC';        
    elseif contains(tmp, 'MLO')
        info.view = 'MLO';        
    else
        info.view = '';
        info.iscc = false;
    end
end
info.view = info.view(~isspace(info.view));
info.ismlo = strcmpi(info.view, 'MLO');
info.iscc = strcmpi(info.view, 'CC');

%SIDE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(full_info,'Laterality')&&~isempty(full_info.Laterality)
    info.side = full_info.Laterality;
elseif isfield(full_info, 'ImageLaterality')
    info.side = full_info.ImageLaterality;
else
    info.side = '';
end

if isfield(full_info, 'Manufacturer')
    info.source = full_info.Manufacturer;
else
    info.source = 'Unknown';
end

if isfield(full_info, 'BodyPartThickness')
    info.H = full_info.BodyPartThickness;
else
    info.H = nan;
end

if isfield(full_info, 'KVP')
    info.KVP = full_info.KVP;
else
    info.KVP = nan;
end

%SPATIAL RESOLUTION (PSIZE)
if isfield(full_info, 'SpatialResolution')
    info.psize = full_info.SpatialResolution;
elseif isfield(full_info, 'DetectorElementPhysicalSize')
    info.psize = full_info.DetectorElementPhysicalSize(1);
elseif isfield(full_info, 'DetectorElementSpacing')
    info.psize = full_info.DetectorElementSpacing(1);
elseif isfield(full_info, 'PixelSpacing')
    info.psize = full_info.PixelSpacing(1);
else
    info.psize = nan;
end

% TARGET MATERIAL
if isfield(full_info, 'AnodeTargetMaterial')
    info.target = full_info.AnodeTargetMaterial(1:2);
else
    info.target = 'NA';
end

% FILTER MATERIAL
if isfield(full_info, 'FilterMaterial')&&~isempty(full_info.FilterMaterial)    
    info.filter = full_info.FilterMaterial(1:2);
else
    info.filter = 'NA';
end

if isfield(full_info, 'StudyDate')
    info.date = full_info.StudyDate;
end

info.iscc = ~info.ismlo;

%COMPRESSION FORCE
if isfield(full_info, 'CompressionForce')&&~isempty(full_info.CompressionForce)
    info.cforce = full_info.CompressionForce;
end

%PATIENT AGE
if isfield(full_info, 'PatientAge')&&~isempty(full_info.PatientAge)
    info.age = full_info.PatientAge;
else
    info.age = '055Y';
end

%EXPOSURE TIME
if isfield(full_info, 'ExposureTime')&&~isempty(full_info.ExposureTime)
    info.exposure = full_info.ExposureTime;
else
    info.exposure = nan;
end

%RADIATION DOSE
if isfield(full_info, 'OrganDose')&&~isempty(full_info.OrganDose)
    info.dose = full_info.OrganDose;
else
    info.dose = nan;
end

if nargout==0, disp(info); end