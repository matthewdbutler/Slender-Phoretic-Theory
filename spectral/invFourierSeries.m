function [solution,mesh] = invFourierSeries(solution,filament,mesh,npoints)
% Recover the calculated concentration and slip as functions of s and theta
% from the stored solutions in terms of Fourier coefficients

% npoints = desired number of theta-points for solution

if filament.isaxisymmetric == 1 %axisymmetric
    solution.c1 = invFS([0;1;0]*solution.c1,npoints);
    solution.c_tot = invFS(solution.c_tot,npoints);
else %non-axisymmetric
    solution.c0 = invFS(solution.c0,npoints);

    % Also invert activity
    mesh.activity_eval = invFS(mesh.activity_eval,npoints);
    nquad = length(mesh.activity_quad(:,1,1));
    activity_quad = zeros(nquad,length(mesh.seval),npoints);
    for i=1:nquad
        % quadrature stored in order (quad,s,FS) rather than (FS,s)
        activity_quad_FS = permute(mesh.activity_quad(i,:,:),[3 2 1]);
        temp_quad = invFS(activity_quad_FS,npoints);  
        activity_quad(i,:,:) = temp_quad';
    end
    mesh.activity_quad = activity_quad;

end

 solution.vlead_x = invFS(solution.vlead_x,npoints);
 solution.vlead_y = invFS(solution.vlead_y,npoints);
 solution.vlead_z = invFS(solution.vlead_z,npoints);

end

function [y,x] = invFS(fy,npoints)
% Converts the vector of Fourier series coeffs for y into the appropriate
% function with npoints in the x-direction
% Cos coeffs are stored as real parts, sin as imaginary

nF = length(fy(:,1))-1; %Highest mode of Fourier series
neval = length(fy(1,:));

x = linspace(0,2*pi,npoints); % x-coord
n = (0:nF)'; % modes
y = zeros(npoints,neval);

for k=1:neval
    cos_terms = real(fy(:,k)).*cos(n.*x);
    sin_terms = imag(fy(:,k)).*sin(n.*x);

    y(:,k) = sum(cos_terms+sin_terms)';
end

end