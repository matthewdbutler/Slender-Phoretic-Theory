function [filament,x] = geometry_examples(centreline,rhochoice,activity)

%% Centreline
% Next function will normalise to arclength = 2
s = linspace(-1,1,1001);
theta = linspace(0,2*pi,1001)';

if  centreline == 1 %Straight rod

    x = [1;0;0]*s;
    filament.isopen = 1;
    filament.planar = 1;

elseif centreline == 2 %planar torus

    phi = linspace(0,2*pi,1001); %parametrise curve
    a = 1;
    x = a*[cos(phi);sin(phi);0*phi];
    filament.isopen = 0;
    filament.planar = 1;

elseif centreline == 3 %U-shape/circular arc

    ang = 1;
    psi = ang*s;
    Px  = spline(s,cos(psi));
    Py  = spline(s,sin(psi));
    x   = [fnval(fnint(Px),s) ; fnval(fnint(Py),s) ; 0*s];
    filament.isopen = 1;
    filament.planar = 1;

elseif centreline == 4 %trefoil

    t = linspace(0,2*pi,1001);
    x1 = sin(t) + 2*sin(2*t);
    x2 = cos(t) - 2*cos(2*t);
    x3 = -sin(3*t);
    x = [x1*cos(4*pi/3)+x2*sin(4*pi/3);-x1*sin(4*pi/3)+x2*cos(4*pi/3);x3];
    filament.isopen = 0;
    filament.planar = 0;

elseif centreline == 5 %trefoil on torus surface

    phi = linspace(0,2*pi,1001);
    x = [(4+cos(3*phi)).*cos(2*phi); (4+cos(3*phi)).*sin(2*phi) ; sin(3*phi)];
    filament.isopen = 0;
    filament.planar = 0;

elseif centreline == 6 %another trefoil on torus surface

    phi = linspace(0,2*pi,1001);
    x = [(4+cos(20*phi)).*cos(2*phi); (4+cos(20*phi)).*sin(2*phi) ; sin(20*phi)];
    filament.isopen = 0;
    filament.planar = 0;

elseif centreline == 7 %mobius strip edges

    phi = linspace(0,2*pi,1001);
    x = [(4+cos(1*phi)).*cos(2*phi); (4+cos(1*phi)).*sin(2*phi) ; sin(1*phi)];
    filament.isopen = 0;
    filament.planar = 0;

elseif centreline == 8 %random 3D curve

    x = rand(3,4);
    filament.isopen = 1;
    filament.planar = 0; %Unless something vanishingly unlikley happens

end

%% Cross sectional profile
filament.param_s = s;
filament.theta = theta;

if rhochoice == 1

    filament.rho = 1 + 0*s; %constant radius

elseif rhochoice == 2

    filament.rho = sqrt(1-s.^2); %prolate spheroidal

elseif rhochoice == 3

    filament.rho = 1 - (2/pi)*(sin(pi*s) + 0.5*sin(2*pi*s));

end

%% Activity
if activity == 1 %Uniform activity

    filament.activity = -1 + 0*s;
    filament.isaxisymmetric = 1;

elseif activity == 2 %Janus activity along centreline

    % Active in s>0 (or -0.5<s<0.5 if looped to avoid end problems)
    filament.activity = 0*s;

    if filament.isopen == 1
        filament.activity(s>0) = -1;
    elseif filament.isopen == 0 % looped - shift so jump is not at \pm 1
        filament.activity(s>-0.5) = -1;
        filament.activity(s>0.5) = 0;
    end

    filament.isaxisymmetric = 1;

elseif activity == 3 % Janus: Upper half active (0<theta<pi) with A=-1

    filament.activity = -heaviside(pi-theta) .* (1+0.*s);
    filament.isaxisymmetric = 0;

elseif activity == 4 % sin(theta) variation

    filament.activity = sin(theta) .* (1+0.*s);
    filament.isaxisymmetric = 0;


elseif activity == 5 %sqrt(1 - s^2)

    filament.activity = sqrt(1-s.^2);
    filament.isaxisymmetric = 1;


elseif activity == 6 % Rotating Janus in theta, activity in <theta<

    filament.activity = -heaviside(cos(pi*s-theta));
    filament.isaxisymmetric = 0;

end

% %% Mobility (potential to add different mobilities at future date)
%
% if mobility == 1 %Uniform mobility
%
filament.mobility = 1 + 0*s;
%
% end

%% Slenderness
filament.epsilon = 0.01;

end