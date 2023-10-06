%Planar uniform torus

%% Geometry
% filament is details on the slender geometry 
[filament,x] = geometry_examples(2,1,1);
filament = generate_geometry_from_xpts(x,filament);

%% Mesh
% mesh is details of surface geometry
nelm  = 20; %Number of elements
nquad = 16;
th_res = 21;
mesh  = generate_mesh(filament,nelm,nquad,th_res,0);

%% Solutions
% Calculate SPT solutions
solution = phoretic_concentration(filament,mesh);

%Slip  and swimming
solution = phoretic_slip(filament,mesh,solution);
solution = stokes_swimming(filament,mesh,solution);

%% Code Validation
% Analytic solution (slip is theta-averaged version)
analytic = analytic_examples(4,1,mesh,filament);

% Compare analytical and numerical SPT calculations
c0_error = 100*(analytic.c0 - solution.c0)./analytic.c0;

vlead = [solution.vlead_x(1,:); solution.vlead_y(1,:); solution.vlead_z(1,:)];
v_avg_error = 100*(analytic.v0 - vecnorm(vlead))./analytic.v0;

solution = invFourierSeries(solution,filament,mesh,th_res);
c1_error = 100*(analytic.c1*ones(1,nelm) - solution.c1)./(analytic.c1*ones(1,nelm));


%% Plot
figure
plot(mesh.seval,c0_error, 'b', mesh.seval,c1_error, 'r')
xlabel('s')
ylabel('conc % error')
legend('c_0','c_1')

figure
plot(mesh.seval,v_avg_error)
xlabel('s')
ylabel('theta-averaged v_0 % error')