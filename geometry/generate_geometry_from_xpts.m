function filament = generate_geometry_from_xpts(x,filament)

%Generate SPT geometry from an inputted ordered set of points in 3D space

%Calculate approximate arclength, and get refined version
s = arclength_from_pts(x);
s_smooth = linspace(0,max(s),1001);

%Spline and get refined coordinates
[~,x] = multispline(x,s,s_smooth);

%Recalculate arclength
s = arclength_from_pts(x);
s_smooth = linspace(0,max(s),1001);

%Recalculate spline and get refined coordinates
[~,x] = multispline(x,s,s_smooth);
s = s_smooth;

%Shift to semilength 2
x = 2*x/max(s);
s = arclength_from_pts(x) - 1;

%Recalculate spline for our version
[P,x] = multispline(x,s,s);

%Tangent
T = deriv_eval(P,s,1);

%Curvature
dt_ds = deriv_eval(P,s,2);
curvature = sqrt(dt_ds(1,:).^2 + dt_ds(2,:).^2 + dt_ds(3,:).^2);

%Normal 
if max(abs(curvature))>1e-6
    N = dt_ds./([1;1;1]*curvature);
else % filament is straight
    N = [0;1;0]*ones(1,1001);
end

%Binormal
B = cross(T,N);

%Torsion
[PB,~] = multispline(B,s,s);
B_s = deriv_eval(PB,s,1);
torsion = B_s./N;

%sort out end points
torsion(:,1) = torsion(:,2);
torsion(:,end) = torsion(:,end-1);

%Take average values of torsion
torsion = mean(torsion);

% Now be a little careful with planar case
if filament.planar == 1 
    torsion = 0*s; 
end

%Theta i
PP = spline(s,torsion);
th_i = ppval(fnint(PP),s);

%Output filament
filament.x = x;
filament.s = s;
filament.tangent = T;
filament.normal = N;
filament.binormal = B;
filament.curvature = curvature;
filament.torsion = torsion;
filament.th_i = th_i;

end