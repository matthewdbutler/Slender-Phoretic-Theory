%Glazed trefoil

%Geometry
[filament,x] = geometry_examples(4,1,3);
filament = generate_geometry_from_xpts(x,filament);

%Mesh
nelm  = 100; %Number of elements
nquad = 16;
th_res = 21;
nF = 21;
mesh  = generate_mesh(filament,nelm,nquad,th_res,nF);

%solution
solution = phoretic_concentration(filament,mesh);

%Slip  and swimming
solution = phoretic_slip(filament,mesh,solution);
solution = stokes_swimming(filament,mesh,solution);

% Recover solution from Fourier series expansion
solution = invFourierSeries(solution,filament,mesh,th_res);


%% Print swimming and rotation speeds
disp('Swimming velocity U =')
disp(solution.U)

disp('Rotational velocity \Omega =')
disp(solution.Omega)