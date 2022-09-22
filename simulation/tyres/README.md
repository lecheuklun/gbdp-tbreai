# Github appendices code index

## Appendix F1. Pacejka-Lite fitting script
B1965run18.mat — TTC lateral data for front tyre, used with TireDataAnalysis.mlx  
B1965run24.mat — TTC lateral data for rear tyre, used with TireDataAnalysis.mlx  
RsqCalc.m — calculates  for Pacejka-Lite  
TireDataAnalysis.mlx — fitting script for Pacejka-Lite  
errorResponseSurf.m — calculates error for response surface variables  
lateral_tire_test.mat — lateral test data for fitting script  
longitudinal_tire_test.mat — longitudinal test data for fitting script  
Functions — helper functions  

## Appendix F2. Pacejka-Lite coefficients
B1965run18_MagicFormula_datapoints.mat — front tyre coefficients, discrete  
B1965run18_MagicFormula_datasurfaces.mat — front tyre coefficients, continuous  
B1965run24_MagicFormula_datapoints.mat — rear tyre coefficients, discrete  
B1965run24_MagicFormula_datasurfaces.mat — rear tyre coefficients, continuous  
lookupCoeffs — explainer script for how to interpret coefficient data  

## Appendix F3. Pacejka-Lite surface plot script
surfPlotIA.m — surface plot script for changing camber, used after TireDataAnalysis.mlx  
surfPlotP.m — surface plot script for changing pressure, used after TireDataAnalysis.mlx  

## Appendix F4. MF6.1 fitting script
convertData.m — converts TTC data for fitting script  
fittingMF61yOnly.m — lateral data fitting script for Magic Formula 6.1  
mfeval.m — Magic Formula 6.1 output  

## Appendix F5. MF6.1 coefficients
coefficientsR18.mat — front tyre coefficients, stored in a struct  
coefficientsR24.mat — rear tyre coefficients, stored in a struct  
MF6.1 front — front tyre coefficients, stored in a Carmaker .tir file  
MF6.1 rear — rear tyre coefficients, stored in a Carmaker .tir file  
