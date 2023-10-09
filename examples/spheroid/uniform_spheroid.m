%% Geometry
% filament is details on the slender geometry 
[filament,x] = geometry_examples(1,2,1);
filament = generate_geometry_from_xpts(x,filament);

%% Mesh
% mesh is details of surface geometry
nelm  = 100; %Number of elements
nquad = 16;
th_res = 20;
mesh  = generate_mesh(filament,nelm,nquad,th_res,0);

%% Solutions
% Calculate SPT solutions
solution = phoretic_concentration(filament,mesh);

% Slip
solution = phoretic_slip(filament,mesh,solution);

% Swimming speed
solution = stokes_swimming(filament,mesh,solution);

%% Validation
% Analytic solution
analytic = analytic_examples(2,1,mesh,filament);

% Error
conc_error = 100.*(analytic.c0-solution.c0)./analytic.c0;
slip_error = 100.*(analytic.v0-solution.vlead_x(1,:))./analytic.v0;

%% Plot
plot(mesh.seval,conc_error); hold on
plot(mesh.seval,slip_error)
legend('% error in conc','% error in slip')