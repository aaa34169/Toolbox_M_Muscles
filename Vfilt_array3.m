% FUNCTION
% Vfilt_array3.m
%__________________________________________________________________________
%
% PURPOSE
% Filtering of vector 
%
% SYNOPSIS
% Vf = Vfilt_array3(V,f,fc)
%
% INPUT
% V (i.e., vector) 
% f (i.e., sampling frequency)
% fc (i.e., cut frequency)
%
% OUTPUT
% Vf (i.e., vector)
%
% DESCRIPTION
% Filtering, along with the 3rd dimension (i.e., all frames, cf. data
% structure in user guide), of the vector components by a 4th order
% Butterworth with special attention when the vector is a column of an
% homogenous matrix
%__________________________________________________________________________
%
% CALLED FUNCTIONS (FROM 3D INVERSE DYNAMICS TOOLBOX) 
% None
% 
% MATLAB VERSION
% Matlab R2012a (without Signal Processing Toolbox)
%__________________________________________________________________________
%
% CHANGELOG
% Created by Rapha�l Dumas
% March 2010
%
% Modified by Raphael Dumas
% September 2012
% Filtering line by line (any number of line)
%__________________________________________________________________________

function Vf = Vfilt_array3(V,f,fc)

% Butterworth
[af,bf] = butter(4,fc./(f/2));

% % No license available for the Signal Processing Toolbox
% if f == 100
%     if fc == 5
%         Butterworth parameters
%         af = [0.0004, 0.0017, 0.0025, 0.0017, 0.0004];
%         bf = [1.0000, -3.1806, 3.8612, -2.1122, 0.4383];
%     elseif fc == 6
%         af = [0.0008, 0.0032, 0.0048, 0.0032, 0.0008];
%         bf = [1.0000, -3.0176, 3.5072, -1.8476, 0.3708];
%     end
% else
%     display('no butterworth parameters available');
%     return
% end
  
% Initialisation
Vf = [];
% Filtering line by line
for c = 1:size(V,1)
    Vc = filtfilt(af,bf,permute(V(c,1,:),[3,1,2]));
    Vf = [Vf;permute(Vc,[3,2,1])];
end


% function y = filtfilt(b,a,x)
% % Subfunction (cut/pasted form the Signal Processing Toolbox)
% 
% len = size(x,1);   % length of input
% b = b(:).';
% a = a(:).';
% nb = length(b);
% na = length(a);
% nfilt = max(nb,na);
% nfact = 3*(nfilt-1);  % length of edge transients
% 
% % set up filter's initial conditions to remove dc offset problems at the
% % beginning and end of the sequence
% if nb < nfilt, b(nfilt)=0; end   % zero-pad if necessary
% if na < nfilt, a(nfilt)=0; end
% % use sparse matrix to solve system of linear equations for initial conditions
% % zi are the steady-state states of the filter b(z)/a(z) in the state-space
% % implementation of the 'filter' command.
% rows = [1:nfilt-1  2:nfilt-1  1:nfilt-2];
% cols = [ones(1,nfilt-1) 2:nfilt-1  2:nfilt-1];
% data = [1+a(2) a(3:nfilt) ones(1,nfilt-2)  -ones(1,nfilt-2)];
% sp = sparse(rows,cols,data);
% zi = sp \ ( b(2:nfilt).' - a(2:nfilt).'*b(1) );
% % non-sparse:
% % zi = ( eye(nfilt-1) - [-a(2:nfilt).' [eye(nfilt-2); zeros(1,nfilt-2)]] ) \ ...
% %      ( b(2:nfilt).' - a(2:nfilt).'*b(1) );
% 
% % Extrapolate beginning and end of data sequence using a "reflection
% % method".  Slopes of original and extrapolated sequences match at
% % the end points.
% % This reduces end effects.
% y = [2*x(1)-x((nfact+1):-1:2);x;2*x(len)-x((len-1):-1:len-nfact)];
% 
% % filter, reverse data, filter again, and reverse data again
% y = filter(b,a,y,zi*y(1));
% y = y(length(y):-1:1);
% y = filter(b,a,y,zi*y(1));
% y = y(length(y):-1:1);
% 
% % remove extrapolated pieces of y
% y([1:nfact len+nfact+(1:nfact)]) = [];

