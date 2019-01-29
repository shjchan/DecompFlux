function [flux2] = fluxcorrect(ST, flux, lb, ub, solver)
%[flux2] = fluxcorrect(ST,flux,lb,ub,eps5)
%To correct the imbalance in flux distribution 'flux' given the
%stoichiometric matrix 'ST', lower bound 'lb', upper bound 'ub' of the flux
%vector and the tolerance of imbalance 'eps5' (default to be 10e-30)

%if nargin < 5
    %eps5 = 10 ^ (-30);
%end
[Nm, Nr] = size(ST);
good = true;
if strcmpi(solver, 'cplex')
    H = 2 * eye(Nr);
    f = -2 * flux;
    z = (flux == 0);
    ub2 = ub;
    ub2(z) = 0;
    options = cplexoptimset;
    options.TolFun = 10^(-9);
    options.TolRLPFun = 10^(-15);
    %options.TolXInteger = 10^(-15);
    %[flux2,fval,~,out]=cplexqp(H,f,[ST;-ST],eps5*ones(2*Nm,1),[],[],lb,ub2,[],options);
    [flux2,fval,~,out]=cplexqp(H, f, [], [], ST, zeros(Nm, 1), lb, ub2, [], options);
    if out.cplexstatus ~= 1
        good = false;
    end
elseif strcmpi(solver, 'cobra')
    QP = struct();
    QP.F = 2 * eye(Nr);
    QP.c = -2 * flux;
    z = (flux == 0);
    QP.ub = ub;
    QP.ub(z) = 0;
    QP.lb = lb;
    QP.A = ST;
    QP.b = zeros(Nm, 1);
    QP.csense = char('E' * ones(1, Nm));
    QP.osense = 1;
    solution = solveCobraQP(QP);
    flux2 = solution.full;
    fval = solution.obj;
    if solution.stat ~= 1
        good = false;
    end
else
    error('Only Cplex or COBRA is supported.')
end
disp('Correct unbalanced flux.');
if ~good
    disp('No solution!');
else
    ib1 = ST * flux;
    ib2 = ST*flux2;
    fprintf('objective value:         %.4e\n', fval + (flux' * flux));
    fprintf('Sum of imbalance before: %.4e\n', sum(abs(ib1)));
    fprintf('Sum of imbalance now:    %.4e\n', sum(abs(ib2)));
    fprintf('No. of nonzeros before:  %d\n', sum(flux ~= 0));
    fprintf('No. of nonzeros now:     %d\n', sum(flux2 ~=0));
end
end

