function [flux2] = fluxcorrect(ST,flux,lb,ub,eps5)
%[flux2] = fluxcorrect(ST,flux,lb,ub,eps5)
%To correct the imbalance in flux distribution 'flux' given the
%stoichiometric matrix 'ST', lower bound 'lb', upper bound 'ub' of the flux
%vector and the tolerance of imbalance 'eps5' (default to be 10e-30)
if nargin<5
    eps5=10^(-30);
end
[Nm,Nr]=size(ST);
H=2*eye(Nr);
f=-2*flux;
z=(flux==0);
ub2=ub;
ub2(z)=0;
options=cplexoptimset;
options.TolFun=10^(-9);
options.TolRLPFun=10^(-15);
options.TolXInteger=10^(-15);
%options.Simplex='on';
%[flux2,fval,~,out]=cplexqp(H,f,[],[],ST,zeros(Nm,1),lb,ub2,[],options);
[flux2,fval,~,out]=cplexqp(H,f,[ST;-ST],eps5*ones(2*Nm,1),[],[],lb,ub2,[],options);
%disp(out);
disp('Correct unbalanced flux.');
if out.cplexstatus~=1
    disp('No solution!');
else
    ib1=ST*flux;
    ib2=ST*flux2;
    disp({'objective value:' (fval+(flux'*flux))});
    disp({'Sum of imbalance before:' sum(abs(ib1));...
        'Sum of imbalance now:' sum(abs(ib2))});
    disp({'No. of nonzeros before:' sum(flux~=0);...
        'No. of nonzeros now:' sum(flux2~=0)});
end
end

