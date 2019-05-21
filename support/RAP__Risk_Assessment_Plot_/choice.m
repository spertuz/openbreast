function [ choice_out ] = choice( varargin )
%  Combines choice vectors so as to enable the selecting out of a specific group of patients.
%  varargin:  vectors of 1s and 0s indicating the preferences
%  There is probably a more efficient way of doing this, but, hey, it works.

n=nargin;
temp=varargin{1};

    for i=2:n
        temp=temp.*varargin{i};
    end

choice_out=temp;

end

