%% Simulation code of manuscript titled "Combined cloud and electricity portfolio optimization for cloud service providers" 
% Contribution: This paper proposes a holistic portfolio optimization framework for a CSP that manages a group of...
% geo-distributed IDCs to make operational decisions in both cloud and electricity markets.
% Author: Caishan Guo
% Version: June 24, 2024
% Environment: Matlab R2023a with Gurobi 10
% All results are saved as a struct called 'Aresults'
% Note: academic use only

clear all;
clc

%% Input parameters
load('parameters.mat');                     

% model modification after revision 
NetworkTransmission = [4,5,6];                % the network transmission delay of relaying requests to IDC 1,2,3 (ms)
etaRobust = [0.98,0.98,0.98];                 % reserve parameter to enhance robustness
Delay = (etaRobust.*(ones(1,IDCnum)*500-NetworkTransmission))/(1000*60*60*10^(-4));    % actual request service delay requirement (convert to an hour)
etaIT = 1./(serrate-1./Delay);                % request-server coefficiency
etaIT = repmat(etaIT',1,T);                   

tic;
%% Decision variables
% Decision variables in Model 1
workloadAllocation=sdpvar(IDCnum,T);  % allocated workloads in each IDC in each time slot
serverNum=sdpvar(IDCnum,T);           % number of switched on servers in each IDC in each time slot
pCH=sdpvar(IDCnum,T);                 % charging power of each BESS in each time slot
pDIS=sdpvar(IDCnum,T);                % discharging power of each BESS in each time slot
pEW=sdpvar(IDCnum,T);                 % IDC's power purchase in wholesale electricity market in each time slot
SOC=sdpvar(IDCnum,T);                 % energy level of each BESS in each time slot
pIDC=sdpvar(IDCnum,T);                % IDC's power consumption in each time slot
pEP=sdpvar(IDCnum,T);                 % IDC's power sold in local energy market

% Decision variables in Model 2
rCO=sdpvar(O,T);
piCO=sdpvar(1,T);
Lagrant_rCO_left=sdpvar(O,T);
Lagrant_rCO_right=sdpvar(O,T);
b_rCO_left=binvar(O,T);
b_rCO_right=binvar(O,T);
ReCO=sdpvar(1,T);

% Decision variables in Model 3
rCFS=sdpvar(1,T);
rCFB=sdpvar(1,T);
rS=sdpvar(S,T);
rB=sdpvar(B,T);
piCF=sdpvar(1,T);
Lagrant_balance=sdpvar(1,T);
Lagrant_rS_left=sdpvar(S,T);
Lagrant_rS_right=sdpvar(S,T);
b_rS_left=binvar(S,T);
b_rS_right=binvar(S,T);
Lagrant_rB_left=sdpvar(B,T);
Lagrant_rB_right=sdpvar(B,T);
b_rB_left=binvar(B,T);
b_rB_right=binvar(B,T);
Lagrant_rCFS_left=sdpvar(1,T);
Lagrant_rCFS_right=sdpvar(1,T);
b_rCFS_left=binvar(1,T);
b_rCFS_right=binvar(1,T);
Lagrant_rCFB_left=sdpvar(1,T);
Lagrant_rCFB_right=sdpvar(1,T);
b_rCFB_left=binvar(1,T);
b_rCFB_right=binvar(1,T);
b_rCF=binvar(1,T);


%% Constraints
cons=[];
% Constraints in Model 1 (Section 3.1.1~3.1.3)
cons=[cons,sum(workloadAllocation)==sum(rCO)+rCFS-rCFB+rA];      
cons=[cons,etaIT.*workloadAllocation<=serverNum];
cons=[cons,0.1*serverNum<=serverNum<=serverMax];  
cons=[cons,workloadAllocation>=0];

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

cons=[cons,pEW+pRES+pDIS==pIDC+pCH+pEP];                                   
cons=[cons,pIDC==etaIDC1.*workloadAllocation+etaIDC2.*serverNum];
cons=[cons,0<=pEP<=0.5*pEW];   % ensure the model is feasible


% KKT constraints and linearization constraints of Model 2 (Section 4.1)

kkt_lower=[];
kkt_lower=[kkt_lower,ones(O,1)*piCO-2*citaCO.*(rCO-rCO0)+Lagrant_rCO_right-Lagrant_rCO_left==0];

M = 10e30;    
kkt_lower=[kkt_lower,0<=Lagrant_rCO_left<=M*b_rCO_left,
                     0<=rCO-rCO_lo<=M*(1-b_rCO_left),
                     0<=Lagrant_rCO_right<=M*b_rCO_right,
                     0<=rCO_up-rCO<=M*(1-b_rCO_right)];

cons=[cons,ReCO<=sum(ones(O,1)*piCO_lo.*rCO+ones(O,1)*piCO.*rCO_up-ones(O,1)*piCO_lo.*rCO_up)];
cons=[cons,ReCO<=sum(ones(O,1)*piCO_up.*rCO+ones(O,1)*piCO.*rCO_lo-ones(O,1)*piCO_up.*rCO_lo)];
cons=[cons,ReCO>=sum(ones(O,1)*piCO_up.*rCO+ones(O,1)*piCO.*rCO_up-ones(O,1)*piCO_up.*rCO_up)];
cons=[cons,ReCO>=sum(ones(O,1)*piCO_lo.*rCO+ones(O,1)*piCO.*rCO_lo-ones(O,1)*piCO_lo.*rCO_lo)];

% KKT constraints and linearization constraints of Model 3 (Section 4.2)
cons=[cons,sum(rS)+rCFS==sum(rB)+rCFB];  

kkt_upper=[];
kkt_upper=[kkt_upper,piS_offer-ones(S,1)*piCF-Lagrant_rS_left+Lagrant_rS_right==0,
                     -piB_bid+ones(B,1)*piCF-Lagrant_rB_left+Lagrant_rB_right==0,
                     pi_bid-piCF-Lagrant_rCFS_left+Lagrant_rCFS_right==0,
                     -pi_bid+piCF-Lagrant_rCFB_left+Lagrant_rCFB_right==0];

M=10e30;     
kkt_upper=[kkt_upper,0<=Lagrant_rS_left<=M*b_rS_left,
                     0<=rS<=M*(1-b_rS_left),
                     0<=Lagrant_rS_right<=M*b_rS_right,
                     0<=rS_up-rS<=M*(1-b_rS_right),
                     0<=Lagrant_rB_left<=M*b_rB_left,
                     0<=rB<=M*(1-b_rB_left),
                     0<=Lagrant_rB_right<=M*b_rB_right,
                     0<=rB_up-rB<=M*(1-b_rB_right),
                     0<=Lagrant_rCFS_left<=M*b_rCFS_left,
                     0<=rCFS<=M*(1-b_rCFS_left),
                     0<=Lagrant_rCFS_right<=M*b_rCFS_right,
                     0<=rCF_up-rCFS<=M*(1-b_rCFS_right),
                     0<=Lagrant_rCFB_left<=M*b_rCFB_left,
                     0<=rCFB<=M*(1-b_rCFB_left),
                     0<=Lagrant_rCFB_right<=M*b_rCFB_right,
                     0<=rCF_up-rCFB<=M*(1-b_rCFB_right),
                     rCFS<=rCF_up.*b_rCF,
                     rCFB<=rCF_up.*(1-b_rCF)]; 
                     
cons=[cons,kkt_lower,kkt_upper];

%% Objective function
ReCF = -sum(piS_offer.*rS+rS_up.*Lagrant_rS_right)+sum(piB_bid.*rB-rB_up.*Lagrant_rB_right);
CostCom = sum(piCM.*workloadAllocation);           % communication cost
CostEW = sum(piEW.*pEW);                           % wholesale energy purchase cost
ProfitEP = sum(piEP.*pEP);                         % local energy market revenue 
CostBESS = sum(costBE.*(pCH+pDIS));                % BESS operation cost
CloudProfit = ReCO+ReCF-CostCom;                   % cloud profit
EnergyProfit = -CostEW+ProfitEP-CostBESS;          % energy profit
obj = sum(CloudProfit)+sum(EnergyProfit);          % total profit

%% Obtain results
options = sdpsettings('verbose',0,'solver', 'gurobi');  
sol = optimize(cons,-obj,options);   

if sol.problem == 0
disp('succcessful solved');
else
disp('error');
yalmiperror(sol.problem)
end
toc;

Aresults.rCO0 = sum(rCO0);       % Fig.4
Aresults.rCO = sum(value(rCO));  % Fig.4
Aresults.piCO = value(piCO);     % Fig.4
Aresults.rCFS = value(rCFS);     % Fig.5
Aresults.rCFB = -value(rCFB);    % Fig.5
Aresults.piCF = value(piCF);     % Fig.5
Aresults.pRES = pRES;            % Fig.6
Aresults.pEW = value(pEW);       % Fig.6
Aresults.pCH = -value(pCH);      % Fig.6
Aresults.pDIS = value(pDIS);     % Fig.6
Aresults.pEP = -value(pEP);      % Fig.6
Aresults.pIDC = value(pIDC);     % Fig.6
Aresults.piEW = 100*piEW;        % Fig.6
Aresults.piEP = 100*piEP;        % Fig.6
Aresults.workloadAllocation = value(workloadAllocation);
Aresults.totalR = sum(Aresults.workloadAllocation);                     % Fig.7
Aresults.ReCO = sum(value(ReCO));                                       % Table III
Aresults.ReCF = sum(value(ReCF));                                       % Table III
Aresults.CostCom_index = sum(piCM.*Aresults.workloadAllocation,2);      % Table III
Aresults.CostEW_index = sum(Aresults.piEW.*Aresults.pEW,2)/100;         % Table III
Aresults.ProfitEP_index = sum(-Aresults.piEP.*Aresults.pEP,2)/100;      % Table III
Aresults.CostBESS_index = sum(costBE.*(-Aresults.pCH+Aresults.pDIS),2); % Table III
Aresults.CostCom = sum(value(CostCom));                                 % Table III
Aresults.CostEW = sum(value(CostEW));                                   % Table III
Aresults.ProfitEP = sum(value(ProfitEP));                               % Table III
Aresults.CostBESS = sum(value(CostBESS));                               % Table III
Aresults.CloudProfit = sum(value(CloudProfit));                         % Table III
Aresults.EnergyProfit = value(EnergyProfit);                            % Table III
Aresults.TotalProfit = sum(value(CloudProfit))+sum(value(EnergyProfit)); % Table III and IV
Aresults.totalPurchase=sum(sum(Aresults.pEW));
Aresults.TableIII_1 = [Aresults.CostCom_index,Aresults.CostEW_index,Aresults.ProfitEP_index,Aresults.CostBESS_index];
Aresults.TableIII_2 = [Aresults.ReCO,Aresults.ReCF,Aresults.CostCom,Aresults.CostEW,Aresults.ProfitEP,Aresults.CostBESS,Aresults.CloudProfit,sum(Aresults.EnergyProfit),Aresults.TotalProfit,Aresults.totalPurchase];
Aresults.Fig4 = [Aresults.rCO;Aresults.rCO0;Aresults.piCO];
Aresults.Fig5 = [Aresults.rCFS;Aresults.rCFB;Aresults.piCF];
Aresults.Fig6 = [Aresults.pRES;Aresults.pEW;Aresults.pDIS;Aresults.pCH;Aresults.pEP;Aresults.pIDC;piEW*100;piEP*100];
Aresults.Fig7 = [Aresults.totalR;Aresults.EnergyProfit];
