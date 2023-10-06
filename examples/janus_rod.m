% Cylinder with axisymmetric activity, half active and half inert

%% Geometry
% filament is details on the slender geometry 
[filament,x] = geometry_examples(1,1,2);
filament = generate_geometry_from_xpts(x,filament);

%% Mesh
% mesh is details of surface geometry
nelm  = 100; %Number of elements
nquad = 16;
th_res = 11;
mesh  = generate_mesh(filament,nelm,nquad,th_res,0);

%% Solutions
% Calculate SPT solutions
solution = phoretic_concentration(filament,mesh);

% Slip
solution = phoretic_slip(filament,mesh,solution);

% Swimming speed
solution = stokes_swimming(filament,mesh,solution);

% Evaluate Fourier series
[solution,mesh] = invFourierSeries(solution,filament,mesh,th_res);

%% Print swimming and rotation speeds
disp('Swimming velocity U =')
disp(solution.U)

disp('Rotational velocity \Omega =')
disp(solution.Omega)

%% Validation
% Analytic solution
analytic = analytic_examples(1,2,mesh,filament);

% Error
conc_error = 100.*(analytic.c0-solution.c0)./analytic.c0;
slip_error = 100.*(analytic.v0-solution.vlead_x(1,:))./analytic.v0;

%% Plot
figure(1)
subplot(1,2,1)
plot(mesh.seval,solution.c_tot,'b-',mesh.seval,analytic.c0,'b--')
legend('numeric conc','analytic conc')
subplot(1,2,2)
plot(mesh.seval,solution.vlead_x(1,:),'r-',mesh.seval,analytic.v0,'r--')
legend('numeric slip','analytic slip')

figure(2)
plot(mesh.seval,conc_error)
hold on
plot(mesh.seval,slip_error)

legend('% error in conc','% error in slip')