function analytic = analytic_examples(geom_option,activity_option,mesh,filament)
% Calculates analytic solutions for several given examples

nelm = length(mesh.seval);
th_res = length(mesh.theta_eval);

if geom_option == 1 %cylinder

    if activity_option == 1 %uniform rod

        analytic.c0 = -log(4.*(1-mesh.seval.^2)./filament.epsilon.^2)./2;
        analytic.c1 = zeros(th_res,nelm);
        analytic.v0 = mesh.seval./(1-mesh.seval.^2);

    elseif activity_option == 2 %Janus rod

        analytic.c0 = - heaviside(mesh.seval).*log(4.*(1-mesh.seval).*mesh.seval./filament.epsilon.^2)./2;
        % Add s<0 part
        analytic.c0 = analytic.c0 - heaviside(-mesh.seval).*log((mesh.seval-1)./mesh.seval)./2;

        analytic.c1 = 0;
        analytic.v0 = -(1./abs(mesh.seval)-1./abs(mesh.seval-1))./2;

    end

elseif geom_option == 2 %spheroid

    if activity_option == 1 %uniform spheroid

        analytic.c0 = -sqrt(1-mesh.seval.^2).*( log(16.*(1-mesh.seval.^2)./filament.epsilon^2)./2 - 1 ) ...
            - mesh.seval.*asin(mesh.seval);
        analytic.c1 = 0;

        analytic.v0 = mesh.seval .* ( log(16.*(1-mesh.seval.^2)./filament.epsilon^2) - 2)./sqrt(1-mesh.seval.^2)./2 ...
            - asin(mesh.seval);

    elseif activity_option == 2 %Janus spheroid
        % Run Michelin & Lauga solution

    end

elseif geom_option == 3 %arc

    if activity_option == 1 %uniform arc

        ang = 1; %Vary in geometry examples
        log_bit = log(64/(ang^2*filament.epsilon^2)*tan(ang*(1+mesh.seval)/4).*tan(ang*(1-mesh.seval)/4));
        analytic.c0 = -0.5*log_bit;
        analytic.c1 = 0.5*ang*(cos(mesh.theta_eval*ones(1,nelm))).*(3 - ones(th_res,1)*log_bit);
        analytic.v0 = NaN;

    end


elseif geom_option == 4 %torus

    if activity_option == 1 %uniform donut

        % Analytical results
        analytic.c0 = -log(8/(filament.epsilon*pi));
        analytic.c1 = pi*cos(mesh.theta_eval)*(-log(8/(filament.epsilon*pi)) + 3/2);

        % Slip velocity (theta-averaged)
        analytic.v0 = pi.*(log(8/filament.epsilon./pi)-1.5)./2;

    elseif activity_option == 2 %dunked (Janus) donut

        % Shift analytic result by 0.5 to match example
        s_shift = mesh.seval + 0.5;
        s_shift(s_shift>1) = s_shift(s_shift>1) - 2;
        tan1 = tan(pi.*(1-s_shift)./4);
        tan2 = tan(pi.*abs(s_shift)./4);
        
        % s>0 parts
        analytic.c0 = -heaviside(s_shift).*log(64.*tan1.*tan2./pi.^2./filament.epsilon.^2)./2;
        analytic.c1 = -heaviside(s_shift).* pi.*cos(mesh.theta_eval)./2 ...
            .*(log(64.*tan1.*tan2./pi.^2./filament.epsilon.^2) - 3);
        % Add s<0 part
        analytic.c0 = analytic.c0 - heaviside(-s_shift).*log(tan1./tan2)./2;
        analytic.c1 = analytic.c1 - heaviside(-s_shift).*pi.*cos(mesh.theta_eval)./2 ...
            .*log(tan1./tan2);

        analytic.v0 = NaN;

    end


end


% Set value for examples not above
if exist('analytic','var') == 0

    analytic.c0 = NaN;
    analytic.v0 = NaN;

end


end