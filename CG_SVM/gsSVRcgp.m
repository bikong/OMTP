           %%*************************************************************%%
           %%         GRID SEARCH  FOR SUPPORT VECTOR REGRESSION          %%
           %%*************************************************************%%
                              %%  File_name:gsSVRcgp.m %%
                              %%  Author: Bikong       %%
          


function [bestmse,bestc,bestg,] = gsSVRcg(train_result,train,gs_option)

if nargin == 2
    gs_option = struct('cmin', -10, 'cmax', 10, ...
                'gmin',-5,'gmax', 5,'pmin',0.01,'pmax',1,'v', 5, ...
                'cstep',0.5,'gstep',0.5,'pstep',0.01,'msestep',0.05);
end

[X,Y,Z] = meshgrid(gs_option.cmin:gs_option.cstep:gs_option.cmax,...
                gs_option.gmin:gs_option.gstep:gs_option.gmax,...
                gs_option.pmin:gs_option.gstep:gs_option.pmax);
[m,n,r] = size(X);
cgp = zeros(m,n,r);
eps = 10^(-4);
bestc = 0;
bestg = 0;
bestp = 0;
bestmse = Inf; 
basenum = 2;
for k = 1:r
    for i = 1:m
        for j = 1:n
            cmd = [' -v ',num2str(gs_option.v),' -c ',num2str( basenum^X(i,j,k) ), ...
            ' -g ',num2str( basenum^Y(i,j,k) ),' -p ',num2str( basenum^X(i,j,k) ), ' -s 3 '];
            cgp(i,j,k) = svmtrain(train_result, train, cmd);
            
            if cgp(i,j,k) < bestmse
                bestmse = cgp(i,j,k);
                bestc = basenum^X(i,j,k);
                bestg = basenum^Y(i,j,k);
                bestp = basenum^Z(i,j,k);

            end
        
            % if abs( cgp(i,j)-mse )<=eps && bestc > basenum^X(i,j,k)
            %     mse = cgp(i,j,k);
            %     bestc = basenum^X(i,j,k);
            %     bestg = basenum^Y(i,j,k);
            %     bestp = basenum^Z(i,j,k);
            % end
        
    end
end

