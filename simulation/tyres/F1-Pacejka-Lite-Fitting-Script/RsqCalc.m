% RSQCALC calculates R^2 for Pacejka-Lite. Requires TireDataAnalysis.mlx to
% have been run. 

% Initialise
Rsq = zeros(3,3,5);

% Loop through bins
for i=1:3
    for n=1:3
        for m=1:5

            % Transformed data
            [Ft,St]= MagicOutput(F_bar_fzia{i,n,m},S_bar_fzia{i,n,m},...
                                    Mu_surf_IA_P{i}(IA_binvalues(n),FZ_binvalues(m)),...
                                    FZ_binvalues(m),CS_surf_IA_P{i}(IA_binvalues(n),FZ_binvalues(m)),datamode);
            
            % Fit line
            Sm = S_M{i,n,m}';
            Fm = -F_M{i,n,m}';
            
            % Limits of fit line
            SmMin = min(Sm);
            SmMax = max(Sm);
            
            % Trim data to within bounds of fit line
            Ft = Ft(find(St > SmMin & St < SmMax));
            St = St(St > SmMin & St < SmMax);
            
            % Corresponding points of data on fit line
            yFt = Fm(dsearchn(Sm, St));
            yBarFt = mean(Ft);
            
            % Sum of square residuals
            ssr = sum((Ft-yFt).^2);
            % Sum of squares total
            sst = sum((Ft-yBarFt).^2);
            
            % Coefficient of determination
            Rsq(i,n,m) = 1-ssr/sst;

        end
    end
end

meanRsq = mean(Rsq,'all')

