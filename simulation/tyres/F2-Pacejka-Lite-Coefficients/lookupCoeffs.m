%{

READ BEFORE USING

The coefficients B,C,D,E from this lookup table is for the
non-dimensionalised Magic Formula.

Meaning Slip values have to be non-dimensionalised before input, and the
resulting output Force is non-dimensional too.

When calculating non-dimensional force & slip you will need Mu (coefficient
of friction) and CS (cornering stiffness). These are also saved in these
.mat files and can be looked up in a similar way.

See images in zip folder for clarification, or ask Cheuk.

%}

clc; clear; close all;

load B1965run18_MagicFormula_datapoints.mat
load B1965run18_MagicFormula_datasurfaces.mat

%{

bp3d1 contains discrete Fz (normal load) testing values, i.e. tyre was tested at 222N, 445N, etc.
Indexed with m, i.e. m=2 is 445N

bp3d2 is the same for IA (inclination angle aka camber) in degrees. 
Indexed with n, i.e. n=3 is 4 deg camber

bp3d3 is the same for P (pressure) in psi. 
Indexed with i, i.e. i=3 is 14 psi
%}

bp3d1;bp3d2;bp3d3;


%{
DISCRETE VALUES ONLY

Tensor of coefficients for B, indexed by (m,n,i)
m-by-n matrix, i pages

i.e. B_Lookup3d_map(1,2,3) means the coefficient B at m=1 (FZ=222), n=2
(IA=2), i=3 (14psi) is 0.7032.

Repeat for other coefficients C,D,E and Mu/CS. Note that D=1 due to Magic Formula
constraints.
%}

B_Lookup3d_map;
B_Lookup3d_map(1,2,3)

%{
INTERPOLATION ONLY

If you wish to interpolate between discrete FZ or IA values for a given
pressure (ie pressure cannot be interpolated), use the following command:

B_surf_IA_P{i}(IA_input,FZ_input)

eg
B_surf_IA_P{2}(0.5,250) means the coefficient B at i=2 (P=12 psi), IA=0.5
deg and FZ=250N is 0.6279.

Repeat for other coefficients.
%}

B_surf_IA_P;
B_surf_IA_P{2}(0.5,250)
