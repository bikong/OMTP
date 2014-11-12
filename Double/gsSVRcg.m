           %%*************************************************************%%
           %%         GRID SEARCH  FOR SUPPORT VECTOR REGRESSION          %%
           %%*************************************************************%%
                              %%  File_name:gsSVRcg.m %%
                              %%  Author: Bikong      %%
          


function [mse,bestc,bestg] = gsSVRcg(train_result,train,gs_option)

if nargin == 2
    gs_option = struct('cmin', -10, 'cmax', 10, ...
                'gmin',-5,'gmax', 5,'v', 5, ...
                'cstep',0.5,'gstep',0.5,'msestep',0.05);
end

[X,Y] = meshgrid(gs_option.cmin:gs_option.cstep:gs_option.cmax,...
                gs_option.gmin:gs_option.gstep:gs_option.gmax);
[m,n] = size(X);
cg = zeros(m,n);
eps = 10^(-4);
bestc = 0;
bestg = 0;
mse = Inf;
basenum = 2;
for i = 1:m
    for j = 1:n
        cmd = [' -v ',num2str(gs_option.v),' -c ',num2str( basenum^X(i,j) ), ...
        ' -g ',num2str( basenum^Y(i,j) ),' -s 3 -p 0.1'];
        cg(i,j) = svmtrain(train_result, train, cmd);
        
        if cg(i,j) < mse
            mse = cg(i,j);
            bestc = basenum^X(i,j);
            bestg = basenum^Y(i,j);
        end
        
        % if abs( cg(i,j)-mse )<=eps && bestc > basenum^X(i,j)
        %     mse = cg(i,j);
        %     bestc = basenum^X(i,j);
        %     bestg = basenum^Y(i,j);
        % end
        
    end
end

% to draw the acc with different c & g
[cg,ps] = mapminmax(cg,0,1);
figure;
[C,h] = contour(X,Y,cg,mse:gs_option.msestep:0.7);
clabel(C,h,'FontSize',7,'Color','r');
xlabel('log2c','FontSize',10);
ylabel('log2g','FontSize',10);
firstline = 'Grid Search'; 
secondline = ['Optimized c=',num2str(bestc),' g=',num2str(bestg), ...
    ' Minimum MSE=',num2str(mse)];
title({firstline;secondline},'Fontsize',11);
grid on;

figure;
meshc(X,Y,cg);
axis([gs_option.cmin,gs_option.cmax,gs_option.gmin,gs_option.gmax,0,1]);
xlabel('log2c','FontSize',10);
ylabel('log2g','FontSize',10);
zlabel('MSE','FontSize',10);
firstline = 'Grid Search'; 
secondline = ['Optimized c=',num2str(bestc),' g=',num2str(bestg), ...
    ' Minimum MSE=',num2str(mse)];
title({firstline;secondline},'Fontsize',11);
