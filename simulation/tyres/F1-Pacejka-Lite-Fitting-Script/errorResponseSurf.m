% ERRORRESPONSESURF calculates percentage error of all variables generated
% by ResponseSurf from their original discrete points.

%% Mu

for i=1:3
    for n=1:3
        for m=1:5
            mu_compare(m,n,i) = mu_fzia{i,n,m};
        end
    end
end

mu_error = mean( abs( (Mu_Lookup3d_map - mu_compare) ./ mu_compare ) * 100, 'all')

%% CS

for i=1:3
    for n=1:3
        for m=1:5
            CS_compare(m,n,i) = CS_fzia{i,n,m};
        end
    end
end

CS_error = mean( abs( (CS_Lookup3d_map - CS_compare) ./ CS_compare ) * 100, 'all')

%% B

for i=1:3
    for n=1:3
        for m=1:5
            B_compare(m,n,i) = B_fzia{i,n,m};
        end
    end
end

B_error = mean( abs( (B_Lookup3d_map - B_compare) ./ B_compare ) * 100, 'all')

%% C

for i=1:3
    for n=1:3
        for m=1:5
            C_compare(m,n,i) = C_fzia{i,n,m};
        end
    end
end

C_error = mean( abs( (C_Lookup3d_map - C_compare) ./ C_compare ) * 100, 'all')

%% D

for i=1:3
    for n=1:3
        for m=1:5
            D_compare(m,n,i) = D_fzia{i,n,m};
        end
    end
end

D_error = mean( abs( (D_Lookup3d_map - D_compare) ./ D_compare ) * 100, 'all')

%% E

for i=1:3
    for n=1:3
        for m=1:5
            E_compare(m,n,i) = E_fzia{i,n,m};
        end
    end
end

E_error = mean( abs( (E_Lookup3d_map - E_compare) ./ E_compare ) * 100, 'all')