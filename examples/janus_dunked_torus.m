%Planar dunked torus

%% Geometry
% filament is details on the slender geometry 
[filament,x] = geometry_examples(2,1,2);
filament = generate_geometry_from_xpts(x,filament);

%% Mesh
% mesh is details of surface geometry
nelm  = 100; %Number of elements
nquad = 16;
th_res = 21;
mesh  = generate_mesh(filament,nelm,nquad,th_res,0);

%% Solutions
% Calculate SPT solutions
solution = phoretic_concentration(filament,mesh);

%Slip and swimming
solution = phoretic_slip(filament,mesh,solution);
solution = stokes_swimming(filament,mesh,solution);

% Recover solution from Fourier series expansion
[solution,mesh] = invFourierSeries(solution,filament,mesh,th_res);
                       
%% Code Validation
% Analytic solution
analytic = analytic_examples(4,2,mesh,filament);

% Compare analytical and numerical SPT calculations
c0_error = 100.*(analytic.c0 - solution.c0)./analytic.c0;
c1_error = 100.*(analytic.c1 - solution.c1)./(analytic.c1);

%% Plot
figure
hold on
plot(mesh.seval,solution.c0);
plot(mesh.seval,analytic.c0,'--');

figure
plot(mesh.seval,c0_error, 'b', mesh.seval,c1_error, 'r')