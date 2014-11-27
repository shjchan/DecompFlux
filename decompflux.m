function [EFM,FM,Rshut,Reduce,Corr,depth] = decompflux(CbModel,flux,...
                                            options,param)
%[EFM,FM,Rshut,Reduce,Corr,depth] = decompflux(FD,CT,SM,solver,param)
%Given the stoichiometrix matrix, a flux distribution can be decomposed
%into a set of EFM by this function.
%
%Input:
%CbModel is a structure with the following fields:
%   S:    m by n stoichiometric matrix 
%   rev:  binary vector indicating reversibility in Sre, 
%         1 for reversible and 0 for irreversible
%flux:  n by 1 flux distribution vector
%
%Optional input:
%options is structure with the following optional fields:
%obj (default to use CbModel.obj if not provided): 
%   n by 1 objective coefficient vector
%   +ve entries for maximization
%   -ve entries for minimization (same as COBRA model)
%CT: the type of constraints you want to use, input either '1' or '2' (as
%   number) for two different types of constraints in the subproblem DC:
%   1:  p_j<=M(v_j)(a_j) where M is a big number v_j is the flux value of 
%       the vector 'flux'
%   2:  p_j<=M(delta_j)(a_j) where M is a big number and delta_j is 0 
%       if v_j in 'flux' is zero, otherwise 1.
%   It is the first type if you omit the input.
%SM: level of display
%    0: (no display)
%    1: display when an EFM is found (default)
%    2: display at every node
%    3: show further the detail of each solve
%solver: currently supported solver: 'cobra' or 'cplex', default 'cobra'.
%block (default to be empty): 
%       If you have k flux modes to be excluded from the possible solutions,
%       then you have a n by k matrix in which each column has 1 for reactions 
%       with non-zero fluxes in the flux mode and 0 for other entries.
% eg. in a network with three reactions only,
%     if block =
%           [1 0;...
%            1 1;...
%            0 1];  
%     then the flux mode with non-zero fluxes in reaction 1 and 2, and the
%     flux mode with non-zero fluxes in reaction 2 and 3, are being blocked
%     from being the solution of the optimization problem.
%
%**The following three fields are not recommended to change for good
%solutions, but depending on your flux distribution, they may need to be
%changed.
%lb (optional, default 0 for all reactions if not exist):   
%       n by 1 lower bound vector (default zeros if empty)
%ub (optional, default 1000 for all reactions if not exist):   
%       n by 1 upper bound vector for scaling purpose. 
%       It is recommended to use values around 1000 
%       because otherwise Cplex is likely to have tolerance problems and 
%       solution may not be accurate.
%big (optional, default 10^8 for CT=1 and 10^7 for CT=2 if not exist):   
%       the large number M used in the subproblem.
%       Not recommended to be too high due to tolerance of the integer solver 
%       of Cplex. Also being too low will over-constrain the set of 
%       possible EFMs. 10^7 for CT=2 should be quite universally applicable
%       but for CT=1, because it uses the information of the flux values in
%       the flux distribution,0 it may depend on your flux distribution. 
%       In general, the more different in the order of magnitude between 
%       the maximum flux and the minimum flux in the flux distribution,  
%       the larger the number is required.
%
%You may use this template:
%options = struct('obj',obj,'CT',1);
%
%param: optimization parameter, for 'cobra' only
%
%Output:
%EFM: the set of EFM
%FM: the set of flux modes found during computation
%Rshut: reactions being shut down when computing the n-th flux mode
%Reduce: the sum of fluxes (of the flux distribution to be decomposed) that
%       is subtracted after finding the n-th EFM
%Corr:  2 by N cell, if corrections of flux modes have taken place due to 
%       high level of flux imbalance, the no. of flux modes and the flux
%       modes before correction will be recorded in Corr.
%depth: When an EFM is found, the depth of the stack, ie., 'current node'
%       number and the depth of the stack after removing useless flux modes, 
%       will be recorded in this N by 2 matrix 'depth'.

%% Pre-process
%Check inputs
if ~exist('options','var'), options = struct(); end
FD = struct();
FD.S = CbModel.S;
if isfield(CbModel,'rev'), 
    FD.rev = CbModel.rev; 
else
    FD.rev = zeros(size(CbModel.S,2),1);
end
FD.flux = flux;
if isfield(CbModel,'rxns'), FD.rxns = CbModel.rxns; end
if isfield(options,'obj'),   
    FD.obj = -options.obj;
elseif isfield(CbModel,'obj'),   
    FD.obj = -CbModel.c;
else
    FD.obj = ones(size(CbModel.S,2),1);
end
if isfield(options,'CT'), CT = options.CT; else CT = 1; end
if isfield(options,'SM'), SM = options.SM; else SM = 1; end
if isfield(options,'solver')
    solver = lower(options.solver);
    if strcmpi(solver, 'cobra')
        solverGood = ~isempty(which('solveCobraMILP'));
    elseif strcmpi(solver, 'cplex')
        solverGood = ~isempty(which('cplexmilp'));
    else
        solverGood = false;
    end
end
if ~isfield(options,'solver') || ~solverGood
    if ~isempty(which('solveCobraMILP'))
        solver='cobra';
    elseif ~isempty(which('cplexmilp'))
        solver='cplex';
    else
        error('A supported solver (Cplex or Cobra solver) cannot be identified.')
    end
    if ~SM
        disp(['Solver changed to: ' solver])
    end
end

if ~isfield(options,'block'),   FD.block=[]; else   FD.block = options.block; end
if isfield(options,'lb')
    if numel(options.lb) ~= size(CbModel.S,2)
        error('length of options.lb not matched with the stoichiometric matrix')
    end
    FD.lb = options.lb;
else
    FD.lb = zeros(size(CbModel.S,2),1);
    FD.lb(FD.rev ~= 0) = -1000;
end
if isfield(options,'ub')
    if numel(options.ub) ~= size(CbModel.S,2)
        error('length of options.ub not matched with the stoichiometric matrix')
    end
    FD.ub = options.ub;   
else
    FD.ub = 1000* ones(size(CbModel.S,2),1);
end
if numel(flux) ~= size(CbModel.S,2)
    error('length of flux not matched with the stoichiometric matrix')
end
if numel(FD.obj) ~= size(CbModel.S,2)
    error('length of options.obj not matched with the stoichiometric matrix')
end    
if ~isempty(FD.block) && size(FD.block,1) ~= size(CbModel.S,2)
    error('No. of rows in options.block must be equal to the no. of columns of the stoichiometric matrix')
end   
%Convert into an irreversible model
FD0 = FD;
FD = SreToSir(FD0);

%% Choose sub-functions according to constraint type designated
t=tic;
if CT==1
    if ~isfield(options,'big'), FD.big=10^8; else FD.big = options.big; end
    if exist('param','var')
        [EFM,FM,Rshut,Reduce,Corr,depth] = decompfluxCT1(FD.S,FD.flux,FD.big,...
            FD.lb,FD.ub,FD.obj,FD.block,SM,solver,param);
    else
        [EFM,FM,Rshut,Reduce,Corr,depth] = decompfluxCT1(FD.S,FD.flux,FD.big,...
            FD.lb,FD.ub,FD.obj,FD.block,SM,solver);
    end
else
    if ~isfield(options,'big'), FD.big=10^7; else FD.big = options.big; end
    if exist('param','var')
        [EFM,FM,Rshut,Reduce,Corr,depth] = decompfluxCT2(FD.S,FD.flux,FD.big,...
            FD.lb,FD.ub,FD.obj,FD.block,SM,solver,param);
    else
        [EFM,FM,Rshut,Reduce,Corr,depth] = decompfluxCT2(FD.S,FD.flux,FD.big,...
            FD.lb,FD.ub,FD.obj,FD.block,SM,solver);
    end
end

%% Post-processing
%convert back to reversible space
EFM2 = EFM(FD.ir2re,:);
EFM2(FD0.rev ~= 0,:) = EFM2(FD0.rev ~= 0,:) - EFM(~FD.ir2re,:);
EFM = EFM2;
clear EFM2
FM2 = FM(FD.ir2re,:);
FM2(FD0.rev ~= 0,:) = FM2(FD0.rev ~= 0,:) - FM(~FD.ir2re,:);
FM = FM2;
if ~isempty(Corr)
    for j = 1:size(Corr,1)
        for k = 1:size(Corr,2)
            v = Corr{j,k}(FD.ir2re);
            v(FD0.rev ~= 0) = v(FD0.rev ~= 0) - Corr{j,k}(~FD.ir2re);
            Corr{j,k} = v;
        end
    end
end
disp({'No. of EFM:' size(EFM,2);'No. of corrected flux modes:' size(Corr,2)});
toc(t)
end


