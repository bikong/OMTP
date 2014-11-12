       %%*************************************************************************%%
       %%             GENETIC ALGORITHM FOR SUPPORT VECTOR REGRESSION             %%
       %%*************************************************************************%%
                         %%======File_name:gaSVRcgp.m=======%%
                         %%======Author:Bikong===============%%


function [BestMSE, Bestc, Bestg, Bestp, ga_option] = gaSVRcgp(train_label,train_data,ga_option)

%% Parameter Initialization
if nargin == 2
    ga_option = struct('maxgen',70,'sizepop',50,'ggap',0.9,...
                'cbound',[0,120],'gbound',[0.01,2], 'pbound',[0.01,1],'v',3);
end
%% maxgen:Maximum Evolved Generation, Default 200, Usual Value Range [100,500]
%% sizepop: Maximum population of Swarm, Default 20, Usual Value range [20,100]
%% cbound = [cmin,cmax], Parameter c Value Range, Default (0,100)
%% gbound = [gmin,gmax], Parameter g Value Range, Default (0,100)
%% v:SVM Cross Validation parameter, defacult 5
%% PRECI: the length of individuals
%% NIND: number of individuals(chrom)
MAXGEN = ga_option.maxgen;
NIND = ga_option.sizepop;
NVAR = 3;
PRECI = 20;
GGAP = ga_option.ggap;
trace = zeros(MAXGEN,2);

%%%******************************** REVIEW OF REQUIRED PARAMETERS :) *****************************%%%
%%fieldD=[len;lb;ub;code;scale;lbin;ubin]
%%len:binary code length of each value(gene)
%%lb: the lower boundary of value 
%%ub: the upper boundary of value 
%%code: coding method of each value  1:standard binary coding   0: binary Gray coding 
%%scale: 1:log scale    0:arithmatic scale 
%%lbin&ubin: 1:include the value at the boundary of range    0:not include

%%stategy:standard binary code; arithmatic scale; first gene only cover lower boundary,second include both
FieldID = ...
[rep([PRECI],[1,NVAR]); ...
[ga_option.cbound(1),ga_option.gbound(1),ga_option.pbound(1); ...
ga_option.cbound(2),ga_option.gbound(2),ga_option.pbound(2)]; ...
[1,1,1;0,0,0;0,1,1;1,1,1]];

%%%*******************************REVIEW OF GA FUNCTION  :) **********************************%%%%
%% crtbp(nind,lind,basevec)   creates a binary population of given size and structure
%% nind:number of individuals(chromosomes); lind: number of value(gene) in each chromosome
%% crtbp(nind,lind)  default: each value(gene) is 0 or 1;
Chrom = crtbp(NIND,NVAR*PRECI);

gen = 1;
v = ga_option.v;
BestMSE = inf;
Bestc = 0;
Bestg = 0;
Bestp = 0;

%%
cg = bs2rv(Chrom,FieldID);

for nind = 1:NIND
    cmd = ['-v ',num2str(v),' -c ',num2str(cg(nind,1)),' -g ',num2str(cg(nind,2)), ...
    ' -p ', num2str(cg(nind,3)),' -s 3'];
    ObjV(nind,1) = svmtrain(train_label,train_data,cmd);
end

[BestMSE,I] = min(ObjV);  %% return the BestMSE and corresponding Column
Bestc = cg(I,1);
Bestg = cg(I,2);
Bestp = cg(I,3);


while 1  
    FitnV = ranking(ObjV);
    
    SelCh = select('sus',Chrom,FitnV,GGAP);
    SelCh = recombin('xovsp',SelCh,0.7);
    SelCh = mut(SelCh);
    
    cg = bs2rv(SelCh,FieldID);
    for nind = 1:size(SelCh,1)
        cmd = ['-v ',num2str(v),' -c ',num2str(cg(nind,1)), ...
        ' -g ',num2str(cg(nind,2)),' -p ', num2str(cg(nind,3)),' -s 3'];
        ObjVSel(nind,1) = svmtrain(train_label,train_data,cmd);
    end
    
    [Chrom,ObjV] = reins(Chrom,SelCh,1,1,ObjV,ObjVSel);   
    
    [NewBestCVaccuracy,I] = min(ObjV);
    cg_temp = bs2rv(Chrom,FieldID);
    temp_NewBestCVaccuracy = NewBestCVaccuracy;
    
    if NewBestCVaccuracy < BestMSE
       BestMSE = NewBestCVaccuracy;
       Bestc = cg_temp(I,1);
       Bestg = cg_temp(I,2);
       Bestp = cg_temp(I,3);
    end
    
    if abs( NewBestCVaccuracy-BestMSE ) <= 10^(-3) && cg_temp(I,1) < Bestc
       BestMSE = NewBestCVaccuracy;
       Bestc = cg_temp(I,1);
       Bestg = cg_temp(I,2);
       Bestp = cg_temp(I,3);
    end    
    
    trace(gen,1) = min(ObjV);
    trace(gen,2) = sum(ObjV)/length(ObjV);

    
    % if gen >= MAXGEN/2 && ( temp_NewBestCVaccuracy-BestMSE ) <= 10^(-2)
    %     break;
    % end
    
    if gen == MAXGEN
        break;
    end
    gen = gen + 1;
end



figure;
hold on;
trace = round(trace*10000)/10000;
plot(trace(1:gen,1),'r*-','LineWidth',1);
plot(trace(1:gen,2),'o-','LineWidth',1);
legend('Best','Average');
xlabel('Generation','FontSize',10);
ylabel('Fitness Value','FontSize',10);
grid on;
axis auto;

line1 = 'Genetic Algorithm';
%line2 = ['terminal generation = ',num2str(gen)];
line3 = ['Optimized p=',num2str(Bestp),' Minimum MSE=',num2str(BestMSE)];
line2 = ['Optimized c=',num2str(Bestc),' Optimized g=',num2str(Bestg)];
title({line1;line2;line3},'FontSize',11);
hold off;


