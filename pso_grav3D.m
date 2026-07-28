function [results]=pso_grav3D(funobj,model,data,opciones,opfun)
%*******************************************************************************
% PSO_GRAV3D
%
% Purpose:
%   Particle Swarm Optimization solver for inverse problems.
%   This function implements the PSO family of algorithms, including
%   CC, CP, PC, PP, PR, PSO, RC, RP, and RR schemes.
%
% Syntax:
%   results = pso_grav3D(funobj, model, data, opciones, opfun)
%
% Inputs:
%   funobj    - Handle to the objective function defining the forward problem.
%   model     - Structure containing the model parameters and search-space
%               bounds.
%   data      - Structure containing the observed data.
%   opciones  - Structure containing the inversion and PSO options.
%   opfun     - Structure containing problem-dependent parameters.
%
% Output:
%   results   - Structure containing the inversion results, model histories,
%               misfit evolution, swarm statistics, and final populations.
%
% Authors:
%   GPINV Research Group
%   University of Oviedo, Spain
%
% Copyright (c) 2026 GPINV Research Group, University of Oviedo
%
% This repository version preserves the computational logic used in the
% associated gravity-inversion experiments.
%*******************************************************************************
idisp =1;
% initializing the random number generator
rand('state',sum(100*clock));
%----------------------------
% Recovering PSO parameters
%----------------------------
oppso=opciones.pso; %PSO structure
%
esquema=oppso.esquema; % PSO algorithm
talla=oppso.size;
maxiter=oppso.maxiter;
%
proyection=oppso.proyection;
elitism=oppso.elitism;
ccontrol=oppso.ccontrol;
cdelta=oppso.cdelta;
% kind of seed
seed=opciones.inversion.seed;
deltat_option=opfun.deltat_option;
%
distmed=[];
inter=[];
nmejoras=[];
leader.inercia=[];
leader.aclocal=[];
leader.acglobal=[];
leader.deltat=[];
%Check for extra fields in text header
if isfield(opfun,'extra_header')
    extra_header = opfun.extra_header;
else
    extra_header = '';
end
%Check for minimum dispersion and deltat if disperion is lower
if isfield(opfun,'rdisp_cc')
    rdisp_cc = opfun.rdisp_cc;
else
    rdisp_cc = 5;
end
if isfield(opfun,'deltat_cc')
    deltat_cc = opfun.deltat_cc;
else
    deltat_cc = 1.5;
end
%..........................................................................
%                       Search space
%..........................................................................
% low and upper limits in row format
lowlimit=model.lowlimit;
upperlimit=model.upperlimit;
nparam=length(lowlimit);
%--------------------------------------------------------------------------
% PSO parameters
%--------------------------------------------------------------------------
if isequal(esquema,'PSO')
    disp('EXECUTING PSO');
    inercia=oppso.pso.inertia;
    aclocal=oppso.pso.philocal;
    acglobal=oppso.pso.phiglobal;
elseif isequal(esquema,'CC')
    disp('EXECUTING CC');
    inercia=oppso.cc.inertia;
    aclocal=oppso.cc.philocal;
    acglobal=oppso.cc.phiglobal;
elseif isequal(esquema,'CP')
   disp('EXECUTING CP');
    inercia=oppso.cp.inertia;
    aclocal=oppso.cp.philocal;
    acglobal=oppso.cp.phiglobal;
elseif isequal(esquema,'PP')
    disp('EXECUTING PP');
    inercia=oppso.pp.inertia;
    aclocal=oppso.pp.philocal;
    acglobal=oppso.pp.phiglobal;
elseif isequal(esquema,'RR')
    %     inercia=linspace(2.5,7,100);
    %     aclocal=3.0*(inercia-1.5);
    %     acglobal=aclocal;
    inercia=oppso.rr.inertia;
    aclocal=oppso.rr.philocal;
    acglobal=oppso.rr.phiglobal;
    disp('RR execution')
elseif isequal(esquema,'RC')
    %disp('EXECUTING RC');
    inercia=oppso.rc.inertia;
    aclocal=oppso.rc.philocal;
    acglobal=oppso.rc.phiglobal;
elseif isequal(esquema,'RP')
    %disp('EXECUTING RP');
    inercia=oppso.rp.inertia;
    aclocal=oppso.rp.philocal;
    acglobal=oppso.rp.phiglobal;
elseif isequal(esquema,'PR')
    %disp('EXECUTING PR');
    inercia=oppso.pr.inertia;
    aclocal=oppso.pr.philocal;
    acglobal=oppso.pr.phiglobal;
elseif isequal(esquema,'PC')||isequal(esquema,'RN')
    %disp('EXECUTING PC');
    inercia=oppso.pc.inertia;
    aclocal=oppso.pc.philocal;
    acglobal=oppso.pc.phiglobal;
end

% keeping a copy of the parameters
allinercia=inercia;
allaclocal=aclocal;
allacglobal=acglobal;
%
ninercia=length(inercia);
% adapting the constants to the size of the swarm
if ninercia==1
    inercia=inercia*ones(1,talla);
    aclocal=aclocal*ones(1,talla);
    acglobal=acglobal*ones(1,talla);
elseif talla> ninercia
    nveces=floor(talla/ninercia);
    sobrante=talla-ninercia*floor(talla/ninercia);
    iner=[];
    acl=[];
    acg=[];
    for k=1:nveces
        prow=randperm(ninercia);
        iner=[iner inercia(prow)];
        acl=[acl aclocal(prow)];
        acg=[acg acglobal(prow)];
    end
    if sobrante>0
        prow=randperm(ninercia);
        iner=[iner inercia(prow(1:sobrante))];
        acl=[acl aclocal(prow(1:sobrante))];
        acg=[acg acglobal(prow(1:sobrante))];
    end
    inercia=iner;
    aclocal=acl;
    acglobal=acg;
    %elseif talla < ninercia
    % we change randomly the parameters in each iteration
    % see on line 247
end
%--------------------------------------------------------------------------
% We initialize the swarm and the velocities
%--------------------------------------------------------------------------
if isequal(seed,'given')&&(~isempty(model.initial))
    swarm0=model.initial; % has the same structure as model
    talla=size(swarm0,1);              % initial swarm size
    % we put the prior into the first population
elseif isequal(seed,'random')||isempty(model.initial)
    % Random generation on the search space: each model is written in a row of matrix swarm0
    [swarm0]=initialpop(talla,lowlimit,upperlimit); %in row format
    disp('Initial random population generated');
end
%--------------------------------------------------------------------------
% If a prior model is given, it is cast into the swarm in a random position
%--------------------------------------------------------------------------
if ~isempty(opfun.prior.model)
    index=ceil(talla*rand(1));
    swarm0(index,:)=opfun.prior.model;
end
swarm=swarm0;
%--------------------------------------------------------------------------
% Initial velocities are initialized to zero
%--------------------------------------------------------------------------
vel(1:talla,1:nparam)=zeros(talla,nparam);
%
low=lowlimit(:);
upper=upperlimit(:);
varmin=repmat(low',talla,1); %lowlimit matrix
varmax=repmat(upper',talla,1); % upperlimit matrix
%--------------------------------------------------------------------------
% forward problem routine
%--------------------------------------------------------------------------
fwproblem=func2str(funobj);
problema_directo=['[misfit,swarm]=' fwproblem '(swarm,data,opciones,opfun);'];
%--------------------------------------------------------------------------
% ITERATIONS
%--------------------------------------------------------------------------
mejores_globales =[];
% the number of parameters for mejores_locales is nparam+1 to keep the
% function value
mejores_locales =ones(talla,nparam+1)*NaN;
%
contador=0;
historia=[];
error_hist=[];
models_outmin=[];
models_outmax=[];
error_iter=[];
fittest=[];

% iteration
k=1;
sinmejora=0; %number of iterations without improving
interior=100;
if idisp==1
    fprintf('\n');
    fprintf('%s   Dismk   Int_modl   Iter_num   Best_fit\n',extra_header);
end
%
for iter=1:maxiter
    historia=[historia;swarm];
    %----------------------------------------------------------------------
    % forward solution
    %----------------------------------------------------------------------
    eval(problema_directo);
    % we keep the predicted data for further statistics
    % realmisfit is the prediction misfit
    realmisfit=misfit(:);
    error_hist=[error_hist realmisfit'];
    %----------------------------------------------------------------------
    % local best initialisation
    %----------------------------------------------------------------------
    if isequal (iter,1)
        %------------------------------------------------------------------
        % first iteration
        %------------------------------------------------------------------
        local_bestval=misfit;
        valorf0=misfit;
        [filas,columnas]=size(valorf0);
        if columnas==1
            valorf0=valorf0';
        end
        mejores_locales =[swarm0 valorf0'];
        mod_sinmejora=[1:talla];
        mejoras_iter=[];
        [gbest_val,mejor]=min(misfit);
        gbest=swarm0(mejor,:);
        error_iter(1)=realmisfit(mejor);
    else
        %------------------------------------------------------------------
        % rest of iterations
        %------------------------------------------------------------------
        [mejoras_iter]=find(local_bestval>misfit);
        mod_sinmejora=find(local_bestval<=misfit);
        nmejoras=[nmejoras length(mejoras_iter)];
        % gbest search (minimum problem: minimize the error)
        real_iter_bestval=min(realmisfit);
        [iter_bestval,mejor]=min(misfit);
        % changing the gbest because we improved the misfit
        if (iter_bestval-gbest_val) < 0
            % improving: we change the global optimum
            gbest_val=iter_bestval;
            gbest=swarm(mejor,:);
            mejores_globales =[mejores_globales; [gbest,iter_bestval]];
            % we keep the history of fittest
            fittest=[fittest; gbest]; % best individual
            error_fittest(k)=realmisfit(mejor);
            sinmejora=0;
        else %local
            sinmejora=sinmejora+1;
            mejores_globales =[mejores_globales; [gbest,gbest_val]];
        end
        error_iter(iter)=realmisfit(mejor);
        % those which are improved are automatically actualize
        local_bestval(mejoras_iter)=misfit(mejoras_iter);
    end
    %Swarm's center of gravity
    %This computation can be performed here because actually swarm = swarm0 in
    %the first iteration due to the assign made in line 160
    center_of_gravity = mean(swarm);
    %------------------------------------------------------------------
    % Monitoring the median distance of the swarm to the center of gravity
    % to control swarm collapse
    %------------------------------------------------------------------
%     dismk=median(norm(swarm-repmat(gbest,talla,1)))/norm(gbest);
    mod_disp2 = (swarm-repmat(center_of_gravity,talla,1)).^2;
    dismk = median(sqrt(sum(mod_disp2')))/norm(center_of_gravity);
    if iter >1
        leader.inercia=[leader.inercia, inercia(mejor)];
        leader.aclocal=[leader.aclocal, aclocal(mejor)];
        leader.acglobal=[leader.acglobal, acglobal(mejor)];
        distmed=[distmed;dismk];
    else
        distmed(1)=dismk;
    end
    if idisp==1
        fprintf('%8.2f',dismk/distmed(1)*100);
    end
    %----------------------------------------------------------------------
    % we change the paremeters at each iteration if talla< ninercia
    %----------------------------------------------------------------------
    if talla < ninercia
        prow=randperm(ninercia);
        inercia=allinercia(prow(1:talla));
        aclocal=allaclocal(prow(1:talla));
        acglobal=allacglobal(prow(1:talla));
    end
    % local best positions
    mejores_locales(mejoras_iter,1:end-1)=swarm(mejoras_iter,:);
    mejores_locales(mejoras_iter,end)=local_bestval(mejoras_iter); %their misfit
    [iter_bestval,mejor] = min(misfit);
    % elistism or not on the gbest
    if isequal(elitism,0) %non elitist (more exploration)
        gbest=swarm(mejor,:);  % mejor global en la iteraci�n iter
    end
    %
    gbest_val=iter_bestval;
    bestpos(iter,1:nparam+1) = [gbest,iter_bestval];
    % random coefficients for acceleration (different for each parameter of
    % each individual)
    rannum1 = rand(talla,nparam);
    rannum2 = rand(talla,nparam);
    % local and global accelerations
    aclocal1 = rannum1.*repmat(aclocal(:),1,nparam);
    acglobal2 = rannum2.*repmat(acglobal(:),1,nparam);
    % putting these variables on row format
    relajacion=repmat(inercia(:),1,nparam);
    localbest=mejores_locales(:,1:end-1); %the last column is the misfit
    %----------------------------------------------------------------------
    %  Deltat (pasot) modality
    %----------------------------------------------------------------------
    if isequal (deltat_option,'lime')
        %------------------------------------------------------------------
        % lime and sand algorithm
        %------------------------------------------------------------------
        niter=opfun.niter;

        pasot=opfun.delta1;
        if rem(iter,niter)==0
            pasot=opfun.delta2; %change to pasot2
        end
    else
        rango=opfun.deltatmax-opfun.deltatmin;
        if isequal (deltat_option,'rand')
            %--------------------------------------------------------------
            % random generation for each particle coordinate
            %--------------------------------------------------------------
            pasot=rango*rand(talla,nparam)+ones(talla,nparam)*opfun.deltatmin;
        elseif isequal (deltat_option,'randp')
            %--------------------------------------------------------------
            % random generation by particle
            %--------------------------------------------------------------
            pasot=rango*rand(talla,1)+ones(talla,1)*opfun.deltatmin;
            pasot=repmat(pasot,1,nparam);
        elseif isequal (deltat_option,'randc')
            %--------------------------------------------------------------
            % random generation by coordinates
            %--------------------------------------------------------------
            pasot=rango*rand(1,nparam)+ones(1,nparam)*opfun.deltatmin;
            pasot=repmat(pasot,talla,1);
        end
    end
    %--------------------------------------------------------------
    %  control of swarm collapse
    %--------------------------------------------------------------
    % if distance between particles less than the indicated relative dispersion
    % we increase the pasot parameter to increase exploration (disperse the
    % swarm
    if dismk/distmed(1)*100<rdisp_cc
        pasot=deltat_cc;
        aclocal1=-aclocal1;
        acglobal2=-acglobal2;
    end
    %--------------------------------------------------------------
    %  4-POINTS PSO FAMILY ALGORITHMS
    %--------------------------------------------------------------
    %--------------------------------------------------------------
    % RC-PSO
    %--------------------------------------------------------------
    if isequal(esquema, 'RC')&&(iter==1)
        swarm0=swarm;
        vel0=vel;
        % one iteration of PSO
        termino1=vel.*(1+(relajacion-1).*pasot);
        termino2=aclocal1.*(localbest-swarm).*pasot;
        termino3=acglobal2.*(repmat(gbest,talla,1)-swarm).*pasot;
        vel=termino1+termino2+termino3;
        %--------------------------
        % updating new population
        %--------------------------
        swarm= swarm + vel.*pasot;
    elseif isequal(esquema, 'RC')&&(iter>1)
        peso=1./((1-relajacion).*pasot);
        termino1=2*peso.*(vel0-pasot.*(swarm-swarm0));
        termino2=peso.*aclocal1.*(localbest-swarm).*pasot;
        termino3=peso.*acglobal2.*(repmat(gbest,talla,1)-swarm).*pasot;
        vel=termino1+termino2+termino3;
        swarm0=swarm; % the current iteration swarm before updating
        vel0=vel; % the current velocity updated is kept for next iteration
        swarm= swarm0 + 2*vel.*pasot; % the new swarm
    end
    %--------------------------------------------------------------
    % RP-PSO
    %--------------------------------------------------------------
    if isequal(esquema, 'RP')&&(iter<=2)
        vel0=vel;
        % one iteration of PSO
        termino1=vel.*(1+(relajacion-1).*pasot);
        termino2=aclocal1.*(localbest-swarm).*pasot;
        termino3=acglobal2.*(repmat(gbest,talla,1)-swarm).*pasot;
        vel=termino1+termino2+termino3;
        %--------------------------
        % updating new population
        %--------------------------
        swarm= swarm + vel.*pasot;
    elseif isequal(esquema, 'RP')&&(iter>2)
        peso=1./((1-relajacion).*pasot);
        termino1=(vel-vel0);
        termino2=aclocal1.*(localbest-swarm).*pasot;
        termino3=acglobal2.*(repmat(gbest,talla,1)-swarm).*pasot;
        %
        vel0=vel; % the current velocity before updating is kept for next iteration
        vel=peso.*(termino1+termino2+termino3);
        swarm= swarm + vel.*pasot; % the new swarm
    end
    %--------------------------------------------------------------
    % PR-PSO
    %--------------------------------------------------------------
    if isequal(esquema, 'PR')&&(iter==1)
        vel0=vel;
        swarm0=swarm;
        % one iteration of PSO
        termino1=vel.*(1+(relajacion-1).*pasot);
        termino2=aclocal1.*(localbest-swarm).*pasot;
        termino3=acglobal2.*(repmat(gbest,talla,1)-swarm).*pasot;
        vel=termino1+termino2+termino3;
        %--------------------------
        % updating new population
        %--------------------------
        swarm= swarm + vel.*pasot;
    elseif isequal(esquema, 'PR')&&(iter>1)
        peso=(1-relajacion).*pasot;
        termino1=vel-peso.*vel0;
        termino2=aclocal1.*(localbest-swarm).*pasot;
        termino3=acglobal2.*(repmat(gbest,talla,1)-swarm).*pasot;
        vel0=vel; % the current velocity before updating is kept for next iteration
        vel=termino1+termino2+termino3;
        swarm= swarm + vel.*pasot; % the new swarm
    end
    %--------------------------------------------------------------
    % PC-PSO
    %--------------------------------------------------------------
    if isequal(esquema, 'PC')&&(iter==1)
        swarm0=swarm;
        vel0=vel;
        % one iteration of PSO
        termino1=vel.*(1+(relajacion-1).*pasot);
        termino2=aclocal1.*(localbest-swarm).*pasot;
        termino3=acglobal2.*(repmat(gbest,talla,1)-swarm).*pasot;
        vel=termino1+termino2+termino3;
        %--------------------------
        % updating new population
        %--------------------------
        swarm= swarm + vel.*pasot;
    elseif isequal(esquema, 'PC')&&(iter>1)
        peso=(relajacion-1).*pasot;
        termino1=peso.*vel/2+(swarm-swarm0);
        termino2=aclocal1.*(localbest-swarm).*pasot/2;
        termino3=acglobal2.*(repmat(gbest,talla,1)-swarm).*pasot/2;
        %
        vel0=vel; % the current velocity before updating is kept for next iteration
        vel=termino1+termino2+termino3;
        swarm0=swarm;
        swarm= swarm0 + 2*vel.*pasot; % the new swarm
    end
    %--------------------------
    % PSO FAMILIY MEMBERS
    %--------------------------
    % CC uses two different leader iterations
    if isequal(esquema, 'CC')&&(iter>1)
        % second part on velnew that is calculated on the
        % next iteration with new acloal1,aclocal2,etc
        ter2=coef2.*aclocal1.*(localbest-swarm);
        ter3=coef2.*acglobal2.*(repmat(gbest,talla,1)-swarm);
        newpart=newpart+ter2+ter3;
        vel=newpart;
    elseif isequal(esquema, 'CC')&&(iter==1)
        newpart=zeros(talla,nparam); %initialised to zero
    end
    %
    if isequal(esquema, 'PSO')
        termino1=vel.*(1+(relajacion-1).*pasot);
        termino2=aclocal1.*(localbest-swarm).*pasot;
        termino3=acglobal2.*(repmat(gbest,talla,1)-swarm).*pasot;
        velnew=termino1+termino2+termino3;
        vel=velnew;
        %--------------------------
        % updating new population
        %--------------------------
        swarm= swarm + vel.*pasot;
    elseif isequal(esquema, 'CC')
        % first part of newpart which is going to be use in the next
        % iteration. This has to be calculated first to update vel!!
        coef1=(2+(relajacion-1).*pasot)./(2+(1-relajacion).*pasot);
        termino1=vel.*coef1;
        %
        coef2=pasot./(2+(1-relajacion).*pasot);
        termino2=coef2.*aclocal1.*(localbest-swarm);
        %
        termino3=coef2.*acglobal2.*(repmat(gbest,talla,1)-swarm);
        newpart=termino1+termino2+termino3;
        % calculate new swarm
        t1=vel.*(1+0.5*(relajacion-1).*pasot);
        t2=0.5*acglobal2.*(repmat(gbest,talla,1)-swarm).*pasot;
        t3=0.5*aclocal1.*(localbest-swarm).*pasot;
        vel=t1+t2+t3; %actualizacion de swarm
        %--------------------------
        % updating new population
        %--------------------------
        swarm= swarm + vel.*pasot;
    elseif isequal(esquema, 'CP')
        if iter >1
            vel=velnew;
        else
            % initial vels are set to zero
        end
        coef1=1-(aclocal1+acglobal2).*pasot.^2;
        coef2=(1+(1-relajacion).*pasot);
        termino1=(coef1./coef2).*vel;
        coef1=aclocal1.*(localbest-swarm).*pasot;
        termino2=coef1./coef2;
        coef1=acglobal2.*(repmat(gbest,talla,1)-swarm).*pasot;
        termino3=coef1./coef2;
        velnew=termino1+termino2+termino3;
        %--------------------------
        % updating new population
        %--------------------------
        swarm= swarm + vel.*pasot;
    elseif isequal(esquema, 'PP')
        if iter >1
            vel=velnew;
        else
            % initial vels are set to zero
        end
        termino1=vel.*(1+(relajacion-1).*pasot);
        termino2=aclocal1.*(localbest-swarm).*pasot;
        termino3=acglobal2.*(repmat(gbest,talla,1)-swarm).*pasot;
        velnew=termino1+termino2+termino3;
        %--------------------------
        % updating new population
        %--------------------------
        swarm= swarm + vel.*pasot;
    elseif isequal(esquema, 'RR')
        coef=1./(1+(1-relajacion).*pasot+(aclocal1+acglobal2).*pasot.^2);
        termino1=vel.*coef;
        termino2=(aclocal1.*coef).*(localbest-swarm).*pasot;
        termino3=(acglobal2.*coef).*(repmat(gbest,talla,1)-swarm).*pasot;
        vel=termino1+termino2+termino3;
        %--------------------------
        % updating new population
        %--------------------------
        swarm= swarm + vel.*pasot;
    elseif isequal(esquema, 'RN')
        swarm=initialpop(talla,lowlimit,upperlimit);
    end

    % if the prior has to be considered I substitute one particle in the
    % swarm (randomly chosen) by the prior
    if (~isempty(opfun.prior.model))&&(rem(iter,opfun.prior.niter)==0)
        index=ceil(talla*rand(1));
        swarm(index,:)=opfun.prior.model;
    end
    %-------------------------------
    % Proyection over search space
    %-------------------------------
    vmin_away = swarm <= varmin;
    vmin_keep      = swarm >  varmin;
    vmax_away = swarm >= varmax;
    vmax_keep      = swarm <  varmax;
    %----------------------------------------------------------------------
    % counting interior models to know how we are sampling the interior of the
    % search space
    %----------------------------------------------------------------------
    parinterior=(swarm >= varmin) & (swarm <= varmax);
    aux=sum(parinterior');
    nvarin=length(find(aux==nparam));
    interior=nvarin/talla*100;
    %%fprintf('Percentage of interior models %f \n',interior);
    if idisp==1
        fprintf('%11.2f',interior);
    end
    inter(iter)=interior;
    %
    % exterior particles: deltat damping
    %
    parexterior=(swarm <= varmin) | (swarm >= varmax);
    if isequal(ccontrol,1)
        pasot=parexterior*cdelta+(1-parexterior).*pasot;
    end
    %--------------------------------------
    %   Proyection method
    %--------------------------------------
    if isequal(proyection, 'near')
        swarm = (vmin_away.*varmin ) + (vmin_keep.*swarm );
        swarm = (vmax_away.*varmax ) + (vmax_keep.*swarm );
    elseif isequal(proyection, 'far')
        swarm = (vmin_away.*varmax ) + (vmin_keep.*swarm );
        swarm = (vmax_away.*varmin ) + (vmax_keep.*swarm );
    elseif isequal(proyection, 'bounce')
        swarm = (vmin_away.*varmin ) + (vmin_keep.*swarm );
        swarm = (vmax_away.*varmax ) + (vmax_keep.*swarm );
        vel = (vel.*vmin_keep) + (-vel.*vmin_away);
        vel = (vel.* vmax_keep) + (-vel.* vmax_away);
    else
        % no proyection
    end
    outmin=length(vmin_away);
    outmax=length(vmax_away);
    models_outmin=[models_outmin,outmin];
    models_outmax=[models_outmax,outmax];
    % Print result to screen
    if idisp==1
        fprintf('    %3d/%3d%11.5f\n',iter,maxiter,iter_bestval);
    end
end

results.historia=historia;
results.fittest=fittest;
results.localbest=localbest;
% important: parent and parent0 are written in the format they had been
% generated (for rebooting purposes)
results.parent=swarm;
results.parent0=swarm0;
% first and last swarm are stored in the format they had been geneerated as
% opciones.inversion.modellog indicates
results.misfit=misfit;
results.error_hist=error_hist;
results.error_iter=error_iter;
results.stat.nmejoras=nmejoras;
results.stat.models_outmin=models_outmin;
results.stat.models_outmax=models_outmax;
results.stat.models_interior=inter;
results.stat.distmed=distmed;
results.stat.leader=leader;  %%at they had been
% generated (for rebooting purposes)
results.parent=swarm;
results.parent0=swarm;