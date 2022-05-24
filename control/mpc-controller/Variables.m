function Model_params = Variables()

Model_params.a = -30:0.1:30;

Model_params.Lr = 0.6885;
Model_params.Lf = 0.8415;
Model_params.m = 210;
Model_params.Iz = 100;

Model_params.Ffz = (Model_params.Lf/(Model_params.Lr+Model_params.Lf)*...
    Model_params.m*9.81)/1000;

Model_params.Cf = 1.9;
Model_params.Df = Model_params.Ffz*(0*Model_params.Ffz+1100);

Model_params.BCDf = 1100*sin(atan(Model_params.Ffz/10)*2);
Model_params.Bf= Model_params.BCDf/(Model_params.Cf*Model_params.Df);
% Model_params.Bf = 10;



Model_params.Frz = (Model_params.Lr/(Model_params.Lr+Model_params.Lf)*...
    Model_params.m*9.81)/1000;

Model_params.Cr = 1.9;
Model_params.Dr = Model_params.Frz*(0*Model_params.Frz+1100);

% Model_params.BCDr = 500*sin(atan(Model_params.Frz/10)*2);
% Model_params.Br= Model_params.BCDr/(Model_params.Cr*Model_params.Dr);
Model_params.Br = 10;


Model_params.T = 2;
Model_params.dt = 0.05;




