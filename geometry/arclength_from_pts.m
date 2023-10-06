function s = arclength_from_pts(x)

%calculate ds
dx = x(1,2:end) - x(1,1:end-1);
dy = x(2,2:end) - x(2,1:end-1);
dz = x(3,2:end) - x(3,1:end-1);
ds = sqrt(dx.^2 + dy.^2 + dz.^2);

%Arclength
s = [0,cumsum(ds)];

end