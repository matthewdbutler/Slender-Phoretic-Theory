%% Geometry
% filament is details on the slender geometry 
[filament,x] = geometry_examples(1,2,2);
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

% Recover solution from Fourier series expansion
[solution,mesh] = invFourierSeries(solution,filament,mesh,th_res);

%% Print swimming and rotation speeds
disp('Swimming velocity U =')
disp(solution.U)

disp('Rotational velocity \Omega =')
disp(solution.Omega)

% Compare with analytic solution from Michelin & Lauga