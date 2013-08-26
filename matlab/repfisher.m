clear
setpath

% import data

start_date = 1983.5;
end_date = 2008.5;
[dpi, labprodg, hourspc, espcg, pcepcg] = getdata(start_date,end_date);

dataset = [dpi,labprodg,hourspc,espcg,pcepcg];
dataset = dataset-kron((mean(dataset)),ones(size(dataset,1),1));

nobs = size(dataset,1);


% set parameters for the estimation exercies
numvar=size(dataset,2);   % number of dependent variables
varlag = 4;               % number of lags
varlagmat = varlag*ones(numvar);
nperiods = 100;           % number of periods for the IRFs
                
nreps = 1000;             % number of bootstrap repetitions

startpos = 1;             % initial observation
endpos = size(dataset,1); % final observation

% get OLS estimates
[coefs,coverr,errmat] = estimate(varlagmat,zeros(numvar),ones(numvar,1),dataset,startpos,endpos);
[cofb,const_b]=ols2ar(coefs(1:end,:),1,0);

% impose LR restrictions through Cholesky decomposition.
coverr = errmat'*errmat/(size(errmat,1));
Vmat = zeros(numvar);
for i = 1:varlag
    Vmat = Vmat - cofb(:,:,i);
end
Vmat = Vmat+eye(numvar);
Vmatinv = eye(numvar)/Vmat;

[cholmat,pchol] = chol(Vmatinv*coverr*Vmatinv');
A0inv = Vmat*cholmat';

% construct identified responses to first and second shock
shock1 = [-100;0;0;0;0];
shock2 = [0;100;0;0;0];
history1=mkirf(A0inv*shock1,cofb,nperiods);
history2=mkirf(A0inv*shock2,cofb,nperiods);

% construct confidence intervals
isconstant = 1;
iscontemp = 0;
yw = 0;
[irfs1, irfs2] = confint(coefs,errmat,dataset,varlag,nreps,nperiods,isconstant,iscontemp,yw,shock1,shock2);

%% turn growth rates to levels, and compute output measure
pireps1 = zeros(nreps,nperiods);
pireps2 = zeros(nreps,nperiods);
labprodreps1 = zeros(nreps,nperiods);
labprodreps2 = zeros(nreps,nperiods);
hoursreps1 = zeros(nreps,nperiods);
hoursreps2 = zeros(nreps,nperiods);
productionreps1 = zeros(nreps,nperiods);
productionreps2 = zeros(nreps,nperiods);
esreps1 = zeros(nreps,nperiods);
esreps2 = zeros(nreps,nperiods);
creps1 = zeros(nreps,nperiods);
creps2 = zeros(nreps,nperiods);

for i=1:nreps
    pireps1(i,:) = cumsum(irfs1(i,1,:));
    pireps2(i,:) = cumsum(irfs2(i,1,:));
    
    labprodreps1(i,:) = cumsum(irfs1(i,2,:));
    labprodreps2(i,:) = cumsum(irfs2(i,2,:));
    
    hoursreps1(i,:) = irfs1(i,3,:);
    hoursreps2(i,:) = irfs2(i,3,:);
    
    esreps1(i,:) = cumsum(irfs1(i,4,:));
    esreps2(i,:) = cumsum(irfs2(i,4,:));
    
    productionreps1(i,:) = labprodreps1(i,:)+hoursreps1(i,:);
    productionreps2(i,:) = labprodreps2(i,:)+hoursreps2(i,:);
   
    creps1(i,:) = cumsum(irfs1(i,5,:));
    creps2(i,:) = cumsum(irfs2(i,5,:));
end

pireps1 = sort(pireps1);
pireps2 = sort(pireps2);
labprodreps1 = sort(labprodreps1);
labprodreps2 = sort(labprodreps2);
hoursreps1 = sort(hoursreps1);
hoursreps2 = sort(hoursreps2);
productionreps1 = sort(productionreps1);
productionreps2 = sort(productionreps2);
esreps1 = sort(esreps1);
esreps2 = sort(esreps2);
creps1 = sort(creps1);
creps2 = sort(creps2);

%% Figure 1, responses to an IST shock

figure; 

subplot(3,2,1);
plot(cumsum(history1(1,:)))
hold on 
plot(pireps1(ceil(nreps*.05),:),'r:')
plot(pireps1(floor(nreps*.95),:),'r:')
hold off
title('Price of investment')

subplot(3,2,2)
plot(cumsum(history1(2,:)))
hold on 
plot(labprodreps1(ceil(nreps*.05),:),'r:')
plot(labprodreps1(floor(nreps*.95),:),'r:')
hold off

title('Labor productivity')

subplot(3,2,3)
plot(history1(3,:))
hold on 
plot(hoursreps1(ceil(nreps*.05),:),'r:')
plot(hoursreps1(floor(nreps*.95),:),'r:')
hold off
title('Hours per capita')

subplot(3,2,4)
plot(cumsum(history1(2,:))+history1(3,:))
hold on 
plot(productionreps1(ceil(nreps*.05),:),'r:')
plot(productionreps1(floor(nreps*.95),:),'r:')
hold off
title('Output per capita')

subplot(3,2,5)
plot(cumsum(history1(4,:)))
hold on 
plot(esreps1(ceil(nreps*.05),:),'r:')
plot(esreps1(floor(nreps*.95),:),'r:')
hold off
title('Investment per capita')

subplot(3,2,6)
plot(cumsum(history1(5,:)))
hold on 
plot(creps1(ceil(nreps*.05),:),'r:')
plot(creps1(floor(nreps*.95),:),'r:')
hold off
title('Consumption per capita')


%% Figure 2, responses to an MFP shock
figure; 

subplot(3,2,1);
plot(cumsum(history2(1,:)))
hold on 
plot(pireps2(ceil(nreps*.05),:),'r:')
plot(pireps2(floor(nreps*.95),:),'r:')
hold off
title('Price of investment')

subplot(3,2,2)
plot(cumsum(history2(2,:)))
hold on 
plot(labprodreps2(ceil(nreps*.05),:),'r:')
plot(labprodreps2(floor(nreps*.95),:),'r:')
hold off

title('Labor productivity')

subplot(3,2,3)
plot(history2(3,:))
hold on 
plot(hoursreps2(ceil(nreps*.05),:),'r:')
plot(hoursreps2(floor(nreps*.95),:),'r:')
hold off
title('Hours per capita')

subplot(3,2,4)
plot(cumsum(history2(2,:))+history2(3,:))
hold on 
plot(productionreps2(ceil(nreps*.05),:),'r:')
plot(productionreps2(floor(nreps*.95),:),'r:')
hold off

title('Output per capita')
subplot(3,2,5)
plot(cumsum(history2(4,:)))
hold on 
plot(esreps2(ceil(nreps*.05),:),'r:')
plot(esreps2(floor(nreps*.95),:),'r:')
hold off
title('Investment per capita')

subplot(3,2,6)
plot(cumsum(history2(5,:)))
hold on 
plot(creps2(ceil(nreps*.05),:),'r:')
plot(creps2(floor(nreps*.95),:),'r:')
hold off
title('Consumption per capita')

%% compute the population variance decompositions

% express the VAR in companion form
[amat bmat]=ar2companion(cofb,A0inv);

% start by computing the unconditional variance in population
gamma0=doublej(amat,bmat*bmat');
% 
nerrs = numvar;
gamma0_partial = zeros(numvar*varlag,numvar*varlag,nerrs);
for i=1:nerrs
     newbmat=0*bmat;
     newbmat(:,i)=bmat(:,i);    
     gamma0_partial(:,:,i)=doublej(amat,newbmat*newbmat');
end
 
omega0 = 2*pi/32;
omega1 = 2*pi/6;
maxiter = 100000;
var_decomp_bc = abs(frequencyplot(omega0,omega1,maxiter,amat,gamma0,gamma0_partial));

% omega0 = 2*pi/200;
% omega1 = 2*pi/6;
% var_decomp_mc = frequencyplot(omega0,omega1,maxiter,newerrlist,amat,gamma0,gamma0_partial);




%% compute the variance decomposition for output and hours worked
% at business cycle frequencies


% exclude dynamics associated with the initial point
% and the constant term

% isolate the dynamics associated with each of the identified shocks

%sample from the errors with replacement
relhvar1=0;
relhvar2=0;
relhvar3=0;
relprodvar1 = 0;
relprodvar2 = 0;
relprodvar3 = 0;

samplelength = 5000;

errpos = round((nobs-varlag-1)*rand(samplelength,1)+ones(samplelength,1));
errmatrep = errmat(errpos,:);

%shock1
A0=inv(A0inv);
ierrmat = A0*errmatrep';
ierrmat1 = 0*ierrmat;
ierrmat1(1,:) = ierrmat(1,:);
errmat1 = A0inv*ierrmat1;
repdata1 = mkymonte_0(dataset,coefs,errmat1');
% compute the output level
hours1 = repdata1(:,3);
production1=cumsum(repdata1(:,1))+hours1;
consumption1 = cumsum(repdata1(:,5));
investment1 = cumsum(repdata1(:,4));
corr1ci = corrcoef(bpass(consumption1,6,32),bpass(investment1,6,32));
corr1cy = corrcoef(bpass(consumption1,6,32),bpass(production1,6,32));


%shock2
% ierrmat2 = 0*ierrmat;
% ierrmat2(2,:) = ierrmat(2,:);
% errmat2 = A0inv*ierrmat2;
% repdata2 = mkymonte_0(dataset,coefs,errmat2');
% hours2 = repdata2(:,3);
% production2=cumsum(repdata2(:,1))+hours2;
% 
% %shock3 and shock4
% ierrmat3 = 0*ierrmat;
% ierrmat3(3,:) = ierrmat(3,:);
% ierrmat3(4,:) = ierrmat(4,:);
% errmat3 = A0inv*ierrmat3;
% repdata3 = mkymonte_0(dataset,coefs,errmat3');
% hours3 = repdata3(:,3);
% production3=cumsum(repdata3(:,1))+hours3;
% 
% % all shocks
% repdata = mkymonte(dataset,coefs,errmatrep);
% hours = repdata(:,3);
% production =cumsum(repdata(:,1))+hours;
% 
% % compute variance decomposition for hours
% relhvar1=relhvar1+var(bpass(hours1,8,32))/var(bpass(hours,8,32));
% relhvar2=relhvar2+var(bpass(hours2,8,32))/var(bpass(hours,8,32));
% relhvar3=relhvar3+var(bpass(hours3,8,32))/var(bpass(hours,8,32));
% 
% % compute variance decompositions for output
% relprodvar1 = relprodvar1+var(bpass(production1,8,32))/var(bpass(production,8,32));
% relprodvar2 = relprodvar2+var(bpass(production2,8,32))/var(bpass(production,8,32));
% relprodvar3 = relprodvar3+var(bpass(production3,8,32))/var(bpass(production,8,32));
% 

% 
% relhvar1=relhvar1/nreps;
% relhvar2=relhvar2/nreps;
% relhvar3=relhvar3/nreps;
% 
% relprodvar1 = relprodvar1/nreps;
% relprodvar2 = relprodvar2/nreps;
% relprodvar3 = relprodvar3/nreps;



