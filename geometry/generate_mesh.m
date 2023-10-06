function mesh = generate_mesh(filament,nelm,nquad,th_res,nFourier)

%Generate mesh from a filament field 

%Evaluation points
s_ends = linspace(-1,1,nelm+1);
mesh.seval = 0.5*(s_ends(2:end) + s_ends(1:end-1));
 
%Quadrature points for each element - allow for uneven distribution
%Split into [s_end-,seval]+[seval,s_end+]
mesh.squad = zeros(2*nquad,nelm);
mesh.quad_weights = zeros(2*nquad,nelm);
for i = 1:nelm
    % interval [seval,s_end+]
    [mesh.squad(1:nquad,i),mesh.quad_weights(1:nquad,i)] = lgwt(nquad,mesh.seval(i),s_ends(i + 1));

    % interval [s_end-,seval]
    [mesh.squad(nquad+1:2*nquad,i),mesh.quad_weights(nquad+1:2*nquad,i)] = lgwt(nquad,s_ends(i),mesh.seval(i));
end

%Quadrature and evaluation points on mesh
[~,mesh.xquad] = multispline(filament.x,filament.s,mesh.squad);
[~,mesh.xeval] = multispline(filament.x,filament.s,mesh.seval);

%Rho (s) on mesh
mesh.rho_quad = singlespline(filament.rho,filament.param_s,mesh.squad);
mesh.rho_eval = singlespline(filament.rho,filament.param_s,mesh.seval);

%Activity on mesh
if filament.isaxisymmetric == 1
    % Activity independent of theta -- write as vector A(s)
    mesh.activity_quad = singlespline(filament.activity,filament.param_s,mesh.squad);
    mesh.activity_eval = singlespline(filament.activity,filament.param_s,mesh.seval);
elseif filament.isaxisymmetric == 0
    % Decompose into Fourier series
    activity_Fourier = FourierSeries(filament.activity,filament.theta,nFourier); % Calculate Fourier Series up to nFourier order
    mesh.activity_quad = singlespline(activity_Fourier,filament.param_s,mesh.squad);
    mesh.activity_eval = singlespline(activity_Fourier,filament.param_s,mesh.seval);
    mesh.activity_quad = permute(mesh.activity_quad,[2 3 1]); % order of indices is now [quad,s,theta Fseries]
end

%Mobility on mesh
mesh.mobility_quad = singlespline(filament.mobility,filament.param_s,mesh.squad);
mesh.mobility_eval = singlespline(filament.mobility,filament.param_s,mesh.seval);

%Curvature of mesh
mesh.curvature_quad = singlespline(filament.curvature,filament.param_s,mesh.squad);
mesh.curvature_eval = singlespline(filament.curvature,filament.param_s,mesh.seval);

%Mesh normal
[~,mesh.normal_quad] = multispline(filament.normal,filament.param_s,mesh.squad);
[~,mesh.normal_eval] = multispline(filament.normal,filament.param_s,mesh.seval);

%Mesh binormal
[~,mesh.binormal_quad] = multispline(filament.binormal,filament.param_s,mesh.squad);
[~,mesh.binormal_eval] = multispline(filament.binormal,filament.param_s,mesh.seval);

% Mesh tangent
[~,mesh.tangent_quad] = multispline(filament.tangent,filament.param_s,mesh.squad);
[~,mesh.tangent_eval] = multispline(filament.tangent,filament.param_s,mesh.seval);

%Evaluation thetas
mesh.theta_eval = linspace(0,2*pi,th_res)';

%Mesh theta_i
mesh.th_i_quad = singlespline(filament.th_i,filament.param_s,mesh.squad);
mesh.th_i_eval = singlespline(filament.th_i,filament.param_s,mesh.seval);

end

function out = singlespline(f,s1,s2)

P.Px = spline(s1,f);
out = ppval(P.Px,s2);

end