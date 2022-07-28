function mask = tri_window(nx,ny,type)
% Given nx and ny, generate a 2D window (triangle or hanning).

if nargin < 3
	type = 1;
end
if type == 1
	mx = [zeros(ceil(nx/4)*2,1);triang(2*floor(nx/4))];
	mx = wshift('1D',mx,ceil(nx/4));
	my = [zeros(ceil(ny/4)*2,1);triang(2*floor(ny/4))];
	my = wshift('1D',my,ceil(ny/4));
	mask = repmat(mx,[1,ny]) .* repmat(my',[nx,1]);
else
	mx = [zeros(ceil(nx/4)*2,1);hann(2*floor(nx/4))];
	mx = wshift('1D',mx,ceil(nx/4));
	my = [zeros(ceil(ny/4)*2,1);hann(2*floor(ny/4))];
	my = wshift('1D',my,ceil(ny/4));
	mask = repmat(mx,[1,ny]) .* repmat(my',[nx,1]);
end

end

