function [ Normals, Disease ] = FindNandD( bm, choice_N, choice_D )
% *********************  FindNandD  ****************************
%   (c) John W Pickering, 2009
%     Christchurch Kidney Research Group
%     University of Otago Christchurch
%     New Zealand
%
%   Last update:  17 July 2012
%
%	Redistribution and use in source and binary forms, with or without 
%   modification, are permitted provided that the following conditions are met:
%
%   * Redistributions of source code must retain the above copyright 
%     notice, this list of conditions and the following disclaimer.
%   * Redistributions in binary form must reproduce the above copyright 
%     notice, this list of conditions and the following disclaimer in 
%     the documentation and/or other materials provided with the distribution
%
% Attribution to John Pickering.  
% *************************************************************************

% PURPOSE
% Find vector of Normals and a vector of diseased from the vector of
% biomarkers, bm, with vectors choice_N, and choice_D picking out the
% normals and diseased


N=logical(choice_N);
D=logical(choice_D);

Norms=bm(N);
Dis=bm(D);

Disease = Dis(~isnan(Dis));
Normals=Norms(~isnan(Norms));

end

