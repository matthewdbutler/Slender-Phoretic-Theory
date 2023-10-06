Code for solving SPT equations numerically 
Works for both open-ended and looped filaments

Solves SPT equations by decomposing into Fourier modes
Before inversion, solution is stored as [i,j]=[Fourier coeffs in theta, arclength]
Fourier coeffs are stored as imaginary numbers with real parts being cos modes, imaginary parts being sin modes

Examples can be found in examples folder, compared to analytic solutions

Issues to avoid:    ensure activity jumps are at least epsilon away from gridpoints
                    do not have activity jumps at end of loop (\pm 1) -- can shift coords around loop to avoid
                    neval should not be too large (suggest <2/epsilon)

Can run pre-set examples using:
geometry_examples(centreline,rhochoice,activity,mobility)

The different options are:
centreline  = 1     straight
            = 2     torus
            = 3     arc/U
            = 4     trefoil
            = 5     Mobius
            = 6     trefoil (alternate)
            = 7     random
rhochoice   = 1     constant (cylinder/rod)
            = 2     prolate spheroid
            = 3     bumps
activity    = 1     uniform
            = 2     Janus (in s)
            = 3     Janus (in theta)
            = 4     sin(theta)
            = 5     sqrt(1-s^2)
mobility    = 1     uniform (only current option)

e.g. uniform cylinder is (1,1,1,1), Janus prolate spheroid is (1,2,2,1), glazed donut is (2,1,3,1)

%%

Alternative: define own slender filament

GEOMETRY
filament.param_s    = input arclength (row vector)
filament.theta      = input angle around cross-section (column vector)
filament.epsilon    = slenderness (<<1)
x                   = centreline geometry in 3D space
filament.rho        = cross-sectional radius

SURFACE PROPERTIES
filament.activity   = activity, given as [theta;s]
filament.mobility   = 1 (no mobility variation allowed yet)

LOGICAL STATEMENTS
filament.isopen         = 1 for open ends, 0 for loops
filament.planar         = 1 for planar filaments, 0 otherwise
filament.isaxisymmetric = 1 if axisymmetric activity, 0 otherwise

%%

Renormalise arclength, calculate tangent/normal/binormal, curvature, etc
filament = generate_geometry_from_xpts(x,filament)

Define the grid for the mesh
nelm        = Number of elements in s (should be even, less than 2/epsilon)
nquad       = (Half) Number of quadrature points per element
nFourier    = Number of terms in Fourier expansion  (not needed for axisymmetric -- can set as 0)
th_res      = Number of elements resolved in theta in results (needs to be split)

Generate the mesh (evaluation and quadrature points)
Quadrature is split on either side of evaluation points
mesh  = generate_mesh(filament,nelm,nquad,th_res,nFourier);

Calculate the SPT surface concentration using quadrature
solution = phoretic_concentration(filament,mesh);

Calculate slip velocity from concentration gradient
solution = phoretic_slip(filament,mesh,solution);

Calculate swimming speed using Koens & Lauga (2017) slender body theory 
solution = stokes_swimming(filament,mesh,solution);

Generate full solution by summing up Fourier series
solution = invFourierSeries(solution,filament,th_res);

%%%
Outputs
mesh = geometry and chemistry defined on evaluation (eval) and quadrature (quad) points
e.g.    mesh.seval          = arclength evaluation points (1 x neval)
        mesh.tangent_eval   = tangent at evaluation points (3 x neval)
        mesh.theta_eval     = theta evaluation points (th_res x 1)
        mesh.xquad          = centreline position at quadrature points (6nquad x neval)
        mesh.activity_quad  = activity at quadrature points (2nquad x neval)
        mesh.quad_weights   = weighting for quadrature integrals (2nquad x neval)

solution = SPT solutions for concentration, slip and swimming
solution.c0         = leading order surface concentration
solution.c1         = next order surface concentration
solution.c_tot      = sum of important concentration contributions, c0+eps*c1
solution.vlead_x    = leading slip contribution in x-direction
solution.vlead_y    = leading slip contribution in y-direction
solution.vlead_z    = leading slip contribution in z-direction
solution.U          = swimiming velocity
solution.Omega      = rotational velocity

Note: before invFourierSeries is run, solution is stored as Fourier modes

