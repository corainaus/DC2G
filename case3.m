%% Simulation code of manuscript titled "Combined cloud and electricity portfolio optimization for cloud service providers" 
% case 3: optimize energy portfolio first, then optimize cloud portfolio

clear all;
clc 

%% Input parameters
load('parameters.mat');                     

tic;

% Decision variables
pCH=sdpvar(IDCnum,T);                 % charging power of each BESS in each time slot
pDIS=sdpvar(IDCnum,T);                % discharging power of each BESS in each time slot
pEW=sdpvar(IDCnum,T);                 % IDC's power purchase in wholesale electricity market in each time slot
SOC=sdpvar(IDCnum,T);                 % energy level of each BESS in each time slot
pIDC=sdpvar(IDCnum,T);                % IDC's power consumption in each time slot
pEP=sdpvar(IDCnum,T);                 % IDC's power sold in local energy market

% Constraints
cons=[];

cons=[cons,SOC(:,1)==SOC0'];
cons=[cons,pCH(:,T)==0,pDIS(:,T)==0];
for k=1:IDCnum
    cons=[cons,0<=pCH(k,:)<=ratedPowerBE(k),0<=pDIS(k,:)<=ratedPowerBE(k)];
    % calculation of SOC
    for t=2:T
        cons=[cons,SOC(k,t)==SOC(k,t-1)+(pCH(k,t-1)*etaCH(k)-pDIS(k,t-1)*etaDIS(k))*deltaT];
    end
    cons=[cons,SOCmin(k)<=SOC(k,:)<=SOCmax(k)];
    cons=[cons,SOC(k,T)==SOC(k,1)];
end

cons=[cons,pEW+pRES+pDIS==pCH+pEP];
cons=[cons,0<=pEW<=100];   % ensure the model is feasible
cons=[cons,0<=pEP<=50];   % ensure the model is feasible

% Objective function
CostEW = sum(piEW.*pEW);                           % wholesale energy purchase cost
ProfitEP = sum(piEP.*pEP);                         % local energy market revenue 
CostBESS = sum(costBE.*(pCH+pDIS));                % BESS operation cost
EnergyProfit = -CostEW+ProfitEP-CostBESS;          % energy profit
obj = sum(EnergyProfit);         

% Obtain results
options = sdpsettings('verbose',0,'solver', 'gurobi');  
sol = optimize(cons,-obj,options);   

if sol.problem == 0
disp('succcessful solved');
else
disp('error');
yalmiperror(sol.problem)
end
toc;


% Aresults.rCO = value(rCO);
% Aresults.piCO = value(piCO);
% Aresults.rCFS = value(rCFS);
% Aresults.rCFB = value(rCFB);
% Aresults.rS = value(rS);
% Aresults.rB = value(rB);
% Aresults.piCF = value(piCF);
Aresults.pRES = pRES;
Aresults.pEW = value(pEW);
Aresults.pCH = value(pCH);
Aresults.pDIS = value(pDIS);
Aresults.pEP = value(pEP);
% Aresults.pIDC = value(pIDC);
% Aresults.workloadAllocation = value(workloadAllocation);
% Aresults.totalR = sum(Aresults.workloadAllocation)/10^4;
% Aresults.serverNum = value(serverNum);
Aresults.SOC = value(SOC);
% Aresults.ReCO = sum(value(ReCO));
% Aresults.ReCF = sum(value(ReCF));
% Aresutls.realReCF=sum(sum(Aresutls.piCF.*(Aresutls.rCFS-Aresutls.rCFB)));
% Aresults.CostCom = sum(value(CostCom));
Aresults.CostEW = sum(value(CostEW));
Aresults.ProfitEP = sum(value(ProfitEP));
Aresults.CostBESS = sum(value(CostBESS));
% Aresults.CloudProfit = sum(value(CloudProfit));
Aresults.EnergyProfit = value(EnergyProfit);
Aresults.TotalProfit = sum(value(EnergyProfit));
% Aresults.CostCom_index = sum(piCM.*Aresults.workloadAllocation,2);
Aresults.CostEW_index = sum(piEW.*Aresults.pEW,2);
Aresults.ProfitEP_index = sum(piEP.*Aresults.pEP,2);                         % 本地售电收益 REP
Aresults.CostBESS_index = sum(costBE.*(Aresults.pCH+Aresults.pDIS),2);                 % 储能运行成本 CBE
Aresults.totalPurchase=sum(sum(Aresults.pEW));
Aresults.TableIII_2 = [Aresults.ProfitEP,Aresults.CostBESS,sum(Aresults.EnergyProfit),Aresults.TotalProfit,Aresults.totalPurchase];
