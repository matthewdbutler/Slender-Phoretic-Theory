function solution = phoretic_concentration(filament,mesh)
% For non-axisymmetric activity, only need to calculate leading order conc
% Activity decomposed as Fourier coefficients
% Zeroth order coefficient requires integral evaluation, others multiples

neval = length(mesh.xeval);
nquad = length(mesh.squad(:,1));
nth = length(mesh.activity_eval(:,1))-1; % Truncation of Fourier Series
xind  =         1:nquad;
yind  =   nquad+1:2*nquad;
zind  = 2*nquad+1:3*nquad;

%Once broken into elements, the order over which you integrate the elements
%doesn't matter, so the integral calculation is the same as unlooped case

leading_order_integral = zeros(1,neval);
activity0_eval = mesh.activity_eval(1,:); % Zeroth Fourier mode for non-axi
activity0_quad = mesh.activity_quad(:,:,1); % Zeroth Fourier mode for non-axi
first_order_integral = zeros(1,neval); % first order integral if axisymm

for i = 1:neval %loop over evaluation points
   
    %Quadrature vectors
    % dx = mesh.xquad(xind,:) - mesh.xeval(1,i);
    % dy = mesh.xquad(yind,:) - mesh.xeval(2,i);
    % dz = mesh.xquad(zind,:) - mesh.xeval(3,i);
    dx =  mesh.xeval(1,i) - mesh.xquad(xind,:);
    dy =  mesh.xeval(2,i) - mesh.xquad(yind,:);
    dz =  mesh.xeval(3,i) - mesh.xquad(zind,:);
    
    %R_0 in equation (distance in integrand)
    r2 = dx.^2 + dy.^2 + dz.^2;
    r = sqrt(r2);
    
    %get 'q'$
    if filament.isopen == 0 %closed filament
        
        q1 = mesh.squad - mesh.seval(i);
        q  = q1;
        q(q1 > 1)  = q(q1 > 1)  - 2;
        q(q1 < -1) = q(q1 < -1) + 2;
        
    else 
    
        %|s - s'| denominator for regularisation
        q = mesh.squad - mesh.seval(i);
        
    end
    
    %Leading order   
    integrand = mesh.rho_quad.*activity0_quad./r - mesh.rho_eval(i).*activity0_eval(i)./abs(q);
    leading_order_integral(i) = sum(sum(integrand.*mesh.quad_weights));
    
    if filament.isaxisymmetric == 1
        % Calculate first order conc c1
        r3 = r.^3;

        % Split integrand into cos and sin modes
        % r0_dot_norm = dx.*mesh.normal_quad(xind,:) + dy.*mesh.normal_quad(yind,:) + dz.*mesh.normal_quad(zind,:);
        % r0_dot_binorm = dx.*mesh.binormal_quad(xind,:) + dy.*mesh.binormal_quad(yind,:) + dz.*mesh.binormal_quad(zind,:);
        r0_dot_norm = dx.*mesh.normal_eval(1,i) + dy.*mesh.normal_eval(2,i) + dz.*mesh.normal_eval(3,i);
        r0_dot_binorm = dx.*mesh.binormal_eval(1,i) + dy.*mesh.binormal_eval(2,i) + dz.*mesh.binormal_eval(3,i);

        integrand_cos = mesh.rho_quad.*mesh.activity_quad.*r0_dot_norm./r3 + mesh.rho_eval(i)*mesh.curvature_eval(i)*mesh.activity_eval(i)./(2*abs(q));
        integrand_sin = mesh.rho_quad.*mesh.activity_quad.*r0_dot_binorm./r3;

        % Combine and integrate (using quadrature weights)
        first_order_integral(i) = -mesh.rho_eval(i)*sum(sum((integrand_cos + 1i*integrand_sin).*mesh.quad_weights));

    end


end

%Leading order log bit
logfactor = 1 - filament.isopen*mesh.seval.^2;
log_leading = mesh.rho_eval.*activity0_eval.*log(4*logfactor./(filament.epsilon*mesh.rho_eval).^2);


if filament.isaxisymmetric == 1
    % Contributions from dc1/dtheta and dc0/ds

    %First-order log bit
    log_first = 0.5*(log(4*logfactor./(filament.epsilon*mesh.rho_eval).^2) - 3).*mesh.rho_eval.^2.*mesh.curvature_eval.*mesh.activity_eval;

    %Output surface concentration
    solution.c0 = 0.5*(leading_order_integral + log_leading);
    solution.c1 = first_order_integral + log_first;

    % Store total solution as first two terms in Fourier series
    % Add in one higher mode (cos2 and sin2) that will be needed later
    solution.c_tot = [solution.c0; filament.epsilon*solution.c1; zeros(1,neval)];

elseif filament.isaxisymmetric == 0
    % Only need to calc theta-deriv of c0 for leading slip

    %Output surface concentration
    solution.c0 = zeros(nth+1,neval);
    solution.c0(1,:)   = 0.5*(leading_order_integral + log_leading);
    solution.c0(2:nth+1,:) = mesh.rho_eval.*mesh.activity_eval(2:nth+1,:)./(1:nth)';

end

end