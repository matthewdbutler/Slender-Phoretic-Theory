function solution = phoretic_slip(filament,mesh,solution)
% Calculate slip velocity from SPT
% theta-dependence stored as cos/sin Fourier coeffs as real/imag part
% e.g. y(n,:) = a+ib represents the a cos(n theta) + b sin(n theta) terms

neval = length(mesh.seval);

% Calculate derivatives in s and theta directions
if filament.isaxisymmetric == 1
    nth = 2;

    %Calculate derivative using spline interpolant
    if filament.isopen == 1 %open ends

        PP = spline(mesh.seval,solution.c0);
        p_der=fnder(PP,1);
        dc0_ds = ppval(p_der,mesh.seval);

    elseif filament.isopen == 0 %looped

        % Account for periodicity by adding an entire copy either side
        c0 = [solution.c0 solution.c0 solution.c0];
        s0 = [mesh.seval-2 mesh.seval mesh.seval+2];

        PP = spline(s0,c0);

        p_der=fnder(PP,1);
        dc0_ds = ppval(p_der,mesh.seval);
    end

    % Should be no theta-dependence in c0
    dc_dth = -1i.*(0:nth)'.*solution.c_tot;

elseif filament.isaxisymmetric == 0
    nth   = length(mesh.activity_eval(:,1))-1;

    % c0 theta-dependent
    % differentiating in theta equivalent to Fourier coeffs x -ni
    dc0_ds = zeros(1,neval); % does not contribute to leading order in slip so ignore
    dc_dth = -1i.*(0:nth)'.*solution.c0;
end

% matrices for multiplying Fourier series by cos or sin
vec = ones(nth,1);
cos_mult_matrix = ( diag(vec,1) + diag(vec,-1) )./2;
sin_mult_matrix = ( -diag(vec,1) + diag(vec,-1) )./2;
cos_mult_matrix(2,1) = 1; sin_mult_matrix(2,1) = 1; % Fix zeroth order contribution

cc = cos_mult_matrix*real(dc_dth); %cos(1)cos(n) result:cos
cs = cos_mult_matrix*imag(dc_dth); %cos(1)sin(n) result:sin
sc = sin_mult_matrix*real(dc_dth); %sin(1)cos(n) result:sin
ss = -sin_mult_matrix*imag(dc_dth); %sin(1)sin(n) result:cos (note sign difference)

% Remove sc and cs first order parts since sincos=sin2 (no const part)
cs(1,:) = 0; sc(1,:) = 0;

%% Calculate leading and first order contributions 
% Sum up cos and sin contributions
vlead_x = (cc+1i*cs).*mesh.binormal_eval(1,:) - (ss+1i*sc).*mesh.normal_eval(1,:);
vlead_y = (cc+1i*cs).*mesh.binormal_eval(2,:) - (ss+1i*sc).*mesh.normal_eval(2,:);
vlead_z = (cc+1i*cs).*mesh.binormal_eval(3,:) - (ss+1i*sc).*mesh.normal_eval(3,:);

% Zeroth order part must be real
vlead_x(1,:) = real(vlead_x(1,:)) + imag(vlead_x(1,:));
vlead_y(1,:) = real(vlead_y(1,:)) + imag(vlead_y(1,:));
vlead_z(1,:) = real(vlead_z(1,:)) + imag(vlead_z(1,:));

solution.vlead_x = vlead_x./mesh.rho_eval./filament.epsilon;
solution.vlead_y = vlead_y./mesh.rho_eval./filament.epsilon;
solution.vlead_z = vlead_z./mesh.rho_eval./filament.epsilon;

% Add on dc0/ds part to zeroth order mode
solution.vlead_x = solution.vlead_x + [1;zeros(nth,1)]*dc0_ds.*mesh.tangent_eval(1,:);
solution.vlead_y = solution.vlead_y + [1;zeros(nth,1)]*dc0_ds.*mesh.tangent_eval(2,:);
solution.vlead_z = solution.vlead_z + [1;zeros(nth,1)]*dc0_ds.*mesh.tangent_eval(3,:);

% Include mobility effect
solution.vlead_x = mesh.mobility_eval.*solution.vlead_x;
solution.vlead_y = mesh.mobility_eval.*solution.vlead_y;
solution.vlead_z = mesh.mobility_eval.*solution.vlead_z;

end
