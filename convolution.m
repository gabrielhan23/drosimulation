function c = convolution(a, b, shape)
%CONVOLUTION MODIFIED BY JULIO CARDENAS, MAY OF 2011.
%   SAME THAN CONV BUT RETURN A VECTOR FROM a(1) to a(end), not the central
%   section as described for the usual convolution function.
%  

if ~isvector(a) || ~isvector(b)
  error(message('MATLAB:conv:AorBNotVector'));
end

if nargin < 3
    shape = 'full';
end

if ~ischar(shape)
  error(message('MATLAB:conv:unknownShapeParameter'));
end

% compute as if both inputs are column vectors
[rows,~]=size(a);
c = conv2(a(:),b(:),shape);
c=c(1:rows);

% restore orientation
if shape(1) == 'f'
    if length(a) > length(b)
        if size(a,1) == 1 %row vector
            c = c.';
        end
    else
        if size(b,1) == 1 %row vector
            c = c.';
        end
    end
else
    if size(a,1) == 1 %row vector
        c = c.';
    end
end

end