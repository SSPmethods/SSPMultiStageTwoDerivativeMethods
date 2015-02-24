%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Optimization Driver File for Finding optimal SSP Two Derivative Runge
%Kutta Methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%===============================================================
%Variable meanings:
% s         - # of stages
% p         - Order of accuracy

% Decision variables:
% A Ahat b bhat - coefficients of the MDRK method
% r             - the radius of absolute monotonicity (multiplied by -1)
% x             - the decision variables stored in the following order:
% x=[A Ahat b bhat r] 
% minreff       -  r  
%===============================================================


%Typically Defined outside m-file, and simply always commented out.
%==============================================
%Starting Types
restart=0;   % Start with randomly generated coeficients
%restart=1;   % Start with loaded method
%restart=2;   % Start with perturbation of loaded method

%Definition of Method
s=2;          %Number of Stages
p=4;          %Order of Method
CC=1/sqrt(2)   %Second Derivative Coefficient (dtVV/dtFE)
minreff =0.01; %Keep looking until method with at least this value is found

%==============================================


%Setting optimization parameters.
%==============================================
if restart==0     
    %Because its generated from random starting point
    %we are being less restrictive on satisfying constraints.
    clear A Ahat b bhat x X r   %Be sure all variables are reset
    X1=[0];
    opts=optimset('MaxFunEvals',100000,'TolCon',1.e-12,'TolFun',1.e-12,'TolX',1.e-12,...
                  'GradObj','on','MaxIter',10000,'Diagnostics','on','Display','on',...
                  'UseParallel','never','Algorithm','sqp');
else

%Restarting from previously found method
%In our finetuning of methods use more restrictive tolerences to be sure 
%our conditions are sharply met
    X1=X;% stores original coefficients from loaded method
    opts=optimset('MaxIter',1000,'MaxFunEvals',10000,...
    'TolCon',1.e-15,'TolFun',1.e-15,'TolX',1.e-15,...
    'GradObj','on','Diagnostics','off','Display','off',...
    'Algorithm','sqp','UseParallel','never');
end

%==============================================



%==============================================
%Now set up:
%The number of unknowns -                     n
%The linear constraints -                     Aeq*x = beq
%The upper and lower bounds on the unknowns - ub, lb

n=(s*(s+1))+1                                %Number of unknowns +1 for R
Aeq=[]; beq=[];                              %No linear equality constraints
bd=10;                                       %Sets non restrictive upper bound on coeficients
lb=zeros(1,n); lb(end)=-2*s;
ub=bd+zeros(1,n); 
ub(end)=0;                                   %requires r>=0
%==============================================
count=0;                                        %Count tracks the number of times optimizer has failed to find a method
info=-2;

%This While loop requires optimizer to keep running until r>minreff is found while satisfying all constraints
while (info==-2 || (r)<minreff || info==0)
    if count==5 %If fails to find a method after 100 times, stop routine
                ('exceed count')
                x=X1;
                r=-x(end);
                break
    end
  %defining initial starting point for Fmincon  
    if restart==1
        x0=X1;          
    elseif restart==2
        x0=X1+.4*rand(1,length(X));
    else 
        x0=[(2*rand(1,n-1)),-.01];         
    end

  %==============================================
    %The optimization call:
    [X,FVAL,info]=fmincon(@mdrk_am_obj,x0,[],[],Aeq,beq,lb,ub,@(x) nlc_mdrk(x,s,p,CC),opts); 
    r=-FVAL;
    count=count+1;
end %while loop
%==============================================
[con,coneq]=nlc_mdrk(X,s,p,CC);
MaxViolation=max([con;coneq']);
[A ,Ahat, b, bhat] =  unpackMSMDRK(X,s);
[Re,P,Q] = Butcher2ShuOsher(A,Ahat,b,bhat,r,CC);
