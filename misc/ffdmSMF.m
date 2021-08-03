function [imn, info] = ffdmSMF(im, info, mask)
% Standardize raw FFDM mammography to SMF [1]
% and normalization)
% Sintax:
%     im = readmammo(dicompath)
%     im = readmammo(dicompath, nflag)
% Inputs:
%     dicompath,  path to dicom file
%     info,       info of mammogram
%     mask,       segmenation mask of breast region
%     
% S. Pertuz
% Nov21/2018
% 
% [1] van Engeland S, Snoeren PR, Huisman H et al (2006) Volumetric
% breast density estimation from full field digital mammo-
% grams. IEEE Trans Med Imaging 25(3):273-W282


% Parameters: %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gap     = 10;       %Gap from breast boundary to inner breast (in mm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if any(size(mask)~=size(im))
    mask = imresize(mask, size(im), 'nearest');
end


% Create table of attenuation coefficients: (Table I in [1])    
[h, KVP] = meshgrid(20:10:90, 24:32);
if strcmpi(info.target, 'RH')&& strcmpi(info.filter, 'RH')
    Du = [.429 .393 .369 .351 .337 .326 .318 .310;...
        .406 .373 .351 .335 .322 .313 .304 .298;...
        .388 .358 .338 .323 .312 .302 .295 .288;...
        .374 .346 .327 .313 .302 .293 .285 .278;...
        .363 .337 .319 .305 .294 .285 .277 .269;...
        .353 .328 .310 .297 .286 .276 .267 .259;...
        .344 .320 .303 .289 .277 .267 .258 .249;...
        .336 .312 .295 .281 .268 .257 .247 .238;...
        .328 .304 .286 .272 .259 .247 .236 .226];
elseif strcmpi(info.target, 'MO')&& strcmpi(info.filter, 'RH')
    Du = [.448 .415 .392 .373 .358 .346 .335 .326;...
        .432 .404 .379 .361 .345 .333 .322 .313;...
        .422 .392 .370 .352 .336 .323 .312 .303;...
        .413 .384 .361 .342 .327 .313 .301 .291;...
        .407 .377 .354 .335 .318 .303 .291 .280;...
        .399 .369 .345 .325 .307 .291 .277 .265;...
        .392 .361 .336 .314 .295 .278 .263 .249;...
        .384 .352 .326 .302 .282 .263 .247 .233;...
        .377 .344 .316 .291 .269 .250 .233 .219];

elseif strcmpi(info.target, 'MO')&& strcmpi(info.filter, 'MO')
    Du = [.513 .477 .452 .433 .417 .403 .390 .379;...
        .497 .761 .435 .414 .396 .379 .363 .348;...
        .484 .448 .421 .398 .378 .358 .340 .324;...
        .468 .432 .403 .377 .354 .333 .313 .295;...
        .456 .419 .388 .360 .335 .312 .291 .273;...
        .442 .403 .370 .340 .313 .289 .267 .249;...
        .429 .387 .352 .320 .292 .267 .246 .228;...
        .414 .371 .334 .301 .272 .247 .227 .210;...
        .402 .357 .319 .285 .256 .232 .212 .196];
else
    error('Data not available for target: %s and filter %s',...
        info.target, info.filter)
end
    
%Evaluate udiff = u_dense - u_fat
if info.KVP>max(KVP(:))
    KVP_ref = max(KVP(:));
elseif info.KVP<min(KVP(:))
    KVP_ref = min(KVP(:));
else
    KVP_ref = info.KVP;
end

if info.H>max(h(:))
    h_ref = max(h(:));
elseif info.H<min(h(:))
    h_ref = min(h(:));
else
    h_ref = info.H;
end

udiff = interp2(h, KVP, Du, h_ref, KVP_ref);

%Obtain reference pixel value of fatty tisse:
mask_inner = mask;
mask_inner(:,1) = false;
mask_inner(1,:) = false;
mask_inner(end,:) = false;
mask_inner(:,end) = false;
mask_inner = imerode(mask_inner, strel('disk', round(0.5*gap/info.psize)));
g_fat = double(max(im(mask_inner)));
    
imn = -log(double(im)/g_fat)/udiff;        
end