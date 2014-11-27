function [FDirr] = SreToSir(FD)
%[Sir,rxnir] = SreToSir(FD)
%
%FD is a structure with the following fields:
%S: stoichiometry matrix with reversible reactions
%rev: binary vector indicating reversibility in Sre, 
%       1 for reversible and 0 for irreversible
%rxns (optional): cell describing names of reactions in S
%flux (optional): flux distribution
%lb (optional): lower bound vector
%ub (optional): upper bound vector
%obj (optional): objective coefficient for Sre
%
%FDirr is a structure with the following fields:
%S: stoichiometry matrix with 
%       each reversible reaction splitted into 2 irreversible reactions
%rxns: cell describing names of reactions in S
%ir2re: a logical vector such that FDirr.S(:,ir2re)=FD.S
%c2: objective coefficient for Sir
rev=FD.rev;
Sre=FD.S;
[Nm,Nr]=size(Sre);
if ~isfield(FD,'flux')
    flux=zeros(Nr,1);
elseif isempty(FD.flux)
    flux=zeros(Nr,1);
else
    flux=FD.flux;
end
if ~isfield(FD,'lb')
    lb=-1000*ones(Nr,1);
elseif isempty(FD.lb)
    lb=-1000*ones(Nr,1);
else
    lb=FD.lb;
end
if ~isfield(FD,'ub')
    ub=1000*ones(Nr,1);
elseif isempty(FD.ub)
    ub=1000*ones(Nr,1);
else
    ub=FD.ub;
end
if ~isfield(FD,'obj')
    c=zeros(Nr,1);
elseif isempty(FD.obj)
    c=zeros(Nr,1);
else
    c=FD.obj;
end
if ~isfield(FD,'rxns')
    FD.rxns=cell(Nr,1);
    for j=1:Nr
        FD.rxns{j}=['R' num2str(j)];
    end
elseif isempty(FD.rxns)
    FD.rxns=cell(Nr,1);
    for j=1:Nr
        FD.rxns{j}=['R' num2str(j)];
    end
end
rev=logical(rev);
Nr2=Nr+sum(rev);
Sro=-Sre(:,rev);
Sir=zeros(Nm,Nr2);
rxnir=cell(Nr2,1);
ir2re=false(Nr2,1);
flux2=zeros(Nr2,1);
lb2=zeros(Nr2,1);
ub2=lb2;
Nrc=0;
Nre=0;
c2=zeros(Nr2,1);
for j=1:Nr
    if rev(j)
        Nrc=Nrc+2;
        Nre=Nre+1;
        Sir(:,(Nrc-1):Nrc)=[Sre(:,j) Sro(:,Nre)];
        rxnir{Nrc-1}=FD.rxns{j};
        rxnir{Nrc}=[FD.rxns{j} '_r'];
        ir2re(Nrc-1)=true;
        flux2(Nrc-1)=max([0 flux(j)]);
        flux2(Nrc)=-min([flux(j) 0]);
        lb2(Nrc-1)=max([0 lb(j)]);
        lb2(Nrc)=-min([ub(j) 0]);
        ub2(Nrc-1)=max([0 ub(j)]);
        ub2(Nrc)=-min([lb(j) 0]);
        c2(Nrc-1:Nrc)=[c(j) -c(j)];
    else
        Nrc=Nrc+1;
        ir2re(Nrc)=true;
        Sir(:,Nrc)=Sre(:,j);
        rxnir{Nrc}=FD.rxns{j};
        flux2(Nrc)=flux(j);
        lb2(Nrc)=lb(j);
        ub2(Nrc)=ub(j);
        c2(Nrc)=c(j);
    end
end
if Nre+Nr~=Nrc || Nre~=sum(rev) || Nrc~=Nr2
    disp({'Nre' Nre; 'Nr' Nr;'Nrc' Nrc;'Nr2' Nr2});
    error('inconsistent dimensions!');
end
if ~isfield(FD,'obj')
   c2=ones(Nr2,1);
elseif isempty(FD.obj)
    c2=ones(Nr2,1);
end
FDirr=FD; %keep the original fields in FD
FDirr.S=Sir;
FDirr.rxns=rxnir;
FDirr.flux=flux2;
FDirr.lb=lb2;
FDirr.ub=ub2;
FDirr.obj=c2;
FDirr.ir2re=ir2re;
end

