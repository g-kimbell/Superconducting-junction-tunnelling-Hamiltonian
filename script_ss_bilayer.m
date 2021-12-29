% Script for investigating properties of superconducting bilayers. Here I
% am investigating the Delta(T) dependence in different superconducting 
% bilayers for different tunnelling parameters:
% weak s-wave / strong s-wave
% weak d-wave / strong d-wave
% weak s-wave / strong d-wave
% I also check if there is any difference between parabolic and
% tight-binding dispersions (there isn't).
%
% Each individual test can take many hours to run. Use the "Run and 
% Advance" (ctrl+shift+enter) or "Run Current Section" (ctrl+enter) to 
% avoid doing every test at once!
%
%% Make new folder for data
mkdir results ss_bilayer

%% Initialise global parameters and define layers

% Global parameters
p=Global_Params();
p.ntest=1000;
p.nfinal=250;
p.abs_tolerance_self_consistency_1S=1e-6;
p.rel_tolerance_Greens=1e-6;

% s-wave high Tc
S1 = Layer;
S1.Delta_0 = 0.0016;
S1.lambda=GKTH_fix_lambda(p,S1,kB*1.764*10);

% s-wave low Tc
S2 = Layer;
S2.Delta_0 = 0.00083;
S2.lambda=GKTH_fix_lambda(p,S2,kB*1.764*5);

% d-wave high Tc
D1 = Layer;
D1.symmetry="d";
D1.Delta_0 = 0.0022;
%D1.lambda=GKTH_fix_lambda(p,D1,0.0023512);
D1.lambda=GKTH_fix_lambda(p,D1,kB*1.764*10*1.32);

% d-wave low Tc
D2 = Layer;
D2.symmetry="d";
D2.Delta_0 = 0.0012;
%D2.lambda=GKTH_fix_lambda(p,D2,0.0010273);
D2.lambda=GKTH_fix_lambda(p,D1,kB*1.764*5*1.32);

% Variables
nTs=50;
Ts=linspace(0.00001,0.001,nTs);
ts=[0,0.25,0.5,1,2,5,10]*1e-3;
nts=length(ts);
Deltas=zeros(2,nts,nTs);

%% s-wave tight-binding dispersion %%
p.nradials=120;
p.lattice_symmetry='mm';
parfor i=1:nts*nTs
    warning('off','all')
    [i1,i2] = ind2sub([nts,nTs],i);
    p1=p;
    p1.ts=ts(i1);
    p1.T=Ts(i2);
    [Ds,~,~] = GKTH_self_consistency_2S_taketurns(p1,[S1,S2]);
    disp("ss: "+i+" of "+nts*nTs)
    Deltas(:,i)=Ds;
end
save("results\ss_bilayer\ss_Delta_T_t_tb.mat",'Ts','Deltas','ts')

%% d-wave tight-binding dispersion %%
p.nradials=60;
p.lattice_symmetry='4mm';
parfor i=1:nts*nTs
    warning('off','all')
    [i1,i2] = ind2sub([nts,nTs],i);
    p1=p;
    p1.ts=ts(i1);
    p1.T=Ts(i2);
    [Ds,~,~] = GKTH_self_consistency_2S_taketurns(p1,[D1,D2]);
    disp("dd: "+i+" of "+nts*nTs)
    Deltas(:,i)=Ds;
end
save("results\ss_bilayer\dd_Delta_T_t_tb.mat",'Ts','Deltas','ts')

% small s-wave/big d-wave tight binding
p.nradials=120;
p.lattice_symmetry='mm';
% S2.phi=pi/2;
% D2.phi=pi/2;
parfor i=1:nts*nTs
    warning('off','all')
    [i1,i2] = ind2sub([nts,nTs],i);
    p1=p;
    p1.ts=ts(i1);
    p1.T=Ts(i2);
    [Ds,~,~] = GKTH_self_consistency_2S_taketurns(p1,[D1,S2]);
    disp("ds: "+i+" of "+nts*nTs)
    Deltas(:,i)=Ds;
end

save("results\ss_bilayer\ds_Delta_T_t_tb.mat",'Ts','Deltas','ts')

% big s-wave/small d-wave tight binding
p.nradials=120;
p.lattice_symmetry='mm';
parfor i=1:nts*nTs
    warning('off','all')
    [i1,i2] = ind2sub([nts,nTs],i);
    p1=p;
    p1.ts=ts(i1);
    p1.T=Ts(i2);
    [Ds,~,~] = GKTH_self_consistency_2S_taketurns(p1,[S1,D2]);
    disp("sd: "+i+" of "+nts*nTs)
    Deltas(:,i)=Ds;
end
save("results\ss_bilayer\sd_Delta_T_t_tb.mat",'Ts','Deltas','ts')

%% small s-wave/big d-wave tight binding
% 
% 
% D1Tc=[0.00094,0.00094,0.00093,0.0009,0.00084,0.00048,0.0004,0.00035];
% S2Tc=[0.00048,0.00048,0.00047,0.00044,0.0003,0.00025,0.000225,0.0002];
% D1D0=[0.0024,0.0024,0.0024,0.0023,0.0023,0.0019,0.0015,0.001];
% S2D0=[0.0008,0.0008,0.0008,0.00075,0.0007,0.0006,0.00055,0.0005];
% 
% D1_guesses=zeros(nts,nTs);
% S2_guesses=zeros(nts,nTs);
% for i = 1:nts
%     D1_guesses(i,:)=D1D0(i).*tanh(1.74.*real(((D1Tc(i)-Ts)./Ts).^0.5))+1e-5;
%     S2_guesses(i,:)=S2D0(i).*tanh(1.74.*real(((S2Tc(i)-Ts)./Ts).^0.5))+1e-5;
% end
% 
% S2.dispersion_type="tb";
% S2.mu=0.06525;
% D1.dispersion_type="tb";
% D1.mu=0.06525;
% D1.lambda=GKTH_fix_lambda(p,D1,0.0023512);
% S2.lambda=GKTH_fix_lambda(p,S2,kB*1.764*5);
% 
% parfor i=1:nts*nTs
%     warning('off','all')
%     [i1,i2] = ind2sub([nts,nTs],i);
%     p1=p;
%     p1.ts=ts(i1);
%     p1.T=Ts(i2);
%     D1_par=D1;
%     S2_par=S1;
%     D1_par.Delta_0=D1_guesses(i1,i2);
%     S2_par.Delta_0=S2_guesses(i1,i2);
%     [Ds,~,~] = GKTH_self_consistency_2S_taketurns(p1,[D1_par,S2_par]);
%     Deltas(:,i)=Ds;
% end
% 
% save("results\ss_bilayer\sd_Delta_T_t_tb_taketurns2.mat",'Ts','Deltas','ts')

%% s-wave parabolic dispersion %%

% S1.dispersion_type="para";
% S2.dispersion_type="para";
% S1.mu=-1.5;
% S2.mu=-1.5;
% 
% S1.lambda=GKTH_fix_lambda(p,S1,kB*1.764*10);
% S2.lambda=GKTH_fix_lambda(p,S2,kB*1.764*5);
% 
% parfor i=1:nts*nTs
%     warning('off','all')
%     [i1,i2] = ind2sub([nts,nTs],i);
%     p1=p;
%     p1.ts=ts(i1);
%     p1.T=Ts(i2);
%     [Ds,~,~] = GKTH_self_consistency_2S(p1,[S1,S2])
%     Deltas(:,i)=Ds;
% end
% save("results\ss_bilayer\ss_Delta_T_t_para.mat",'Ts','Deltas','ts')

%% d-wave parabolic dispersion %%

% D1.dispersion_type="para";
% D2.dispersion_type="para";
% D1.mu=-1.5;
% D2.mu=-1.5;
% 
% D1.lambda=GKTH_fix_lambda(p,D1,0.0038);
% D2.lambda=GKTH_fix_lambda(p,D2,0.0018);
% 
% parfor i=1:nts*nTs
%     warning('off','all')
%     [i1,i2] = ind2sub([nts,nTs],i);
%     p1=p;
%     p1.ts=ts(i1);
%     p1.T=Ts(i2);
%     [Ds,~,~] = GKTH_self_consistency_2S(p1,[D1,D2])
%     Deltas(:,i)=Ds;
% end
% 
% save("results\ss_bilayer\dd_Delta_T_t_para.mat",'Ts','Deltas','ts')

%% small s-wave/big d-wave para

% S2.dispersion_type="para";
% D1.dispersion_type="para";
% S2.mu=-1.5;
% D1.mu=-1.5;
% D1.lambda=GKTH_fix_lambda(p,D1,0.0038);
% S2.lambda=GKTH_fix_lambda(p,S2,kB*1.764*5);
% 
% parfor i=1:nts*nTs
%     warning('off','all')
%     [i1,i2] = ind2sub([nts,nTs],i);
%     p1=p;
%     p1.ts=ts(i1);
%     p1.T=Ts(i2);
%     [Ds,~,~] = GKTH_self_consistency_2S(p1,[D1,S2])
%     Deltas(:,i)=Ds;
% end
% 
% save("results\ss_bilayer\sd_Delta_T_t_para.mat",'Ts','Deltas','ts')

%% big s-wave/small d-wave para

% S1.dispersion_type="para";
% D2.dispersion_type="para";
% S1.mu=-1.5;
% D2.mu=-1.5;
% S1.lambda=GKTH_fix_lambda(p,S1,kB*1.764*10);
% D2.lambda=GKTH_fix_lambda(p,D2,0.0018);
% 
% parfor i=1:nts*nTs
%     warning('off','all')
%     [i1,i2] = ind2sub([nts,nTs],i);
%     p1=p;
%     p1.ts=ts(i1);
%     p1.T=Ts(i2);
%     [Ds,~,~] = GKTH_self_consistency_2S(p1,[S1,D2])
%     Deltas(:,i)=Ds;
% end
% 
% save("results\ss_bilayer\sd_Delta_T_t_para.mat",'Ts','Deltas','ts')

%% Looking at k-space

D2.phi=0;
S2.phi=0;
p.ts=0;
[Fs1,Fk1,kx1,ky1] = GKTH_Greens_radial(p,[D1,S2],verbose=true);
p.ts=2e-3;
[Fs2,Fk2,kx2,ky2] = GKTH_Greens_radial(p,[D1,S2],verbose=true);

figure(1)
subplot(2,2,1)
surf(kx1,ky1,real(squeeze(Fk1(1,:,:))))
shading interp
subplot(2,2,2)
surf(kx1,ky1,real(squeeze(Fk1(2,:,:))))
shading interp
subplot(2,2,3)
surf(kx2,ky2,real(squeeze(Fk2(1,:,:))))
shading interp
subplot(2,2,4)
surf(kx2,ky2,real(squeeze(Fk2(2,:,:))))
shading interp

%%
p=Global_Params();
p.ntest=1000;
p.nfinal=200;
p.nradials=30;
p.lattice_symmetry='4mm';
p.abs_tolerance_self_consistency_1S=1e-6;
p.rel_tolerance_Greens=1e-6;
Deltas=zeros(nTs,1);
D1.lambda=GKTH_fix_lambda(p,D1,kB*1.764*5*1.32);
parfor i=1:nTs
    p1=p;
    p1.T=Ts(i);
    Deltas(i) = GKTH_self_consistency_1S(p1,D1,1,0);
end
hold on
plot(Ts/kB,Deltas,'o-');

%% Eigenvalue plotting
figure(3)
n=1000;
ts=[0,1e-3,1e-2];
p.nkpoints=n;
n=p.nkpoints;
for j=1:length(ts)
    p.ts=ts(j);
    [eigenvalues] = GKTH_find_spectrum(p,[S1,S2]);
    subplot(length(ts),1,j)
    for i=1:8
        plot(linspace(0,1,n),squeeze(eigenvalues(1,:,i)),'.-','color','#0072BD')
        plot(linspace(0,-1,n),diag(eigenvalues(:,:,i)),'.-','color','#0072BD')
        hold on
    end
end

%% Delta_0 calculation
n=24;
ts=[logspace(-5,-3.25,n),logspace(-3.25,-2.75,n),logspace(-2.75,-1,n)];
Deltas=zeros(2,3*n);
p.nfinal=250;
p.nradials=120;
p.lattice_symmetry='mm';
p.T=0.1*kB;
parfor i=1:3*n
    tic
    p1=p;
    p1.ts=ts(i);
    [Ds,~,~] = GKTH_self_consistency_2S_taketurns(p1,[D1,S2]);
    Deltas(:,i)=Ds;
    disp("ds: "+i+" of "+3*n)
    toc
end
save("results\ss_bilayer\ds_Delta0_t.mat",'ts','Deltas')
figure(1)
semilogx(ts,Deltas(1,:),'o-')
hold on
semilogx(ts,Deltas(2,:),'o-')

n=24;
ts=[logspace(-5,-3.25,n),logspace(-3.25,-2.75,n),logspace(-2.75,-1,n)];
Deltas=zeros(2,3*n);
p.nfinal=250;
p.nradials=120;
p.lattice_symmetry='mm';
p.T=0.1*kB;
parfor i=1:3*n
    tic
    p1=p;
    p1.ts=ts(i);
    [Ds,~,~] = GKTH_self_consistency_2S_taketurns(p1,[S1,D2]);
    Deltas(:,i)=Ds;
    disp("sd: "+i+" of "+3*n)
    toc
end
save("results\ss_bilayer\sd_Delta0_t.mat",'ts','Deltas')
figure(2)
semilogx(ts,Deltas(1,:),'o-')
hold on
semilogx(ts,Deltas(2,:),'o-')

n=24;
ts=[logspace(-5,-3.25,n),logspace(-3.25,-2.75,n),logspace(-2.75,-1,n)];
Deltas=zeros(2,3*n);
p.nfinal=250;
p.nradials=120;
p.lattice_symmetry='mm';
p.T=0.1*kB;
parfor i=1:3*n
    tic
    p1=p;
    p1.ts=ts(i);
    [Ds,~,~,~] = GKTH_self_consistency_iterate(p1,[D1,S2]);
    Deltas(:,i)=Ds;
    disp("ds iterate: "+i+" of "+3*n)
    toc
end
save("results\ss_bilayer\ds_Delta0_t_iter.mat",'ts','Deltas')
figure(3)
semilogx(ts,Deltas(1,:),'o-')
hold on
semilogx(ts,Deltas(2,:),'o-')

n=24;
ts=[logspace(-5,-3.25,n),logspace(-3.25,-2.75,n),logspace(-2.75,-1,n)];
Deltas=zeros(2,3*n);
p.nfinal=250;
p.nradials=120;
p.lattice_symmetry='mm';
p.T=0.1*kB;
parfor i=1:3*n
    tic
    p1=p;
    p1.ts=ts(i);
    [Ds,~,~,~] = GKTH_self_consistency_iterate(p1,[S1,D2]);
    Deltas(:,i)=Ds;
    disp("sd iterate: "+i+" of "+3*n)
    toc
end
save("results\ss_bilayer\sd_Delta0_t_iter.mat",'ts','Deltas')
figure(4)
semilogx(ts,Deltas(1,:),'o-')
hold on
semilogx(ts,Deltas(2,:),'o-')

%% Delta_0 calculation
n=48;
ts=linspace(8.25e-4,9.75e-4,n);
Deltas=zeros(2,n);
p.nfinal=250;
p.nradials=120;
p.lattice_symmetry='mm';
p.T=0.1*kB;
parfor i=1:n
    tic
    p1=p;
    p1.ts=ts(i);
    [Ds,~,~] = GKTH_self_consistency_2S_taketurns(p1,[D1,S2]);
    Deltas(:,i)=Ds;
    disp("ds: "+i+" of "+n)
    toc
end
save("results\ss_bilayer\ds_Delta0_t_zoom.mat",'ts','Deltas')
figure(1)
hold on
semilogx(ts,Deltas(1,:),'o-')
semilogx(ts,Deltas(2,:),'o-')

n=48;
ts=linspace(8.25e-4,9.75e-4,n);
Deltas=zeros(2,n);
p.nfinal=250;
p.nradials=120;
p.lattice_symmetry='mm';
p.T=0.1*kB;
parfor i=1:n
    tic
    p1=p;
    p1.ts=ts(i);
    [Ds,~,~] = GKTH_self_consistency_2S_taketurns(p1,[S1,D2]);
    Deltas(:,i)=Ds;
    disp("sd: "+i+" of "+n)
    toc
end
save("results\ss_bilayer\sd_Delta0_t_zoom.mat",'ts','Deltas')
figure(2)
hold on
semilogx(ts,Deltas(1,:),'o-')
semilogx(ts,Deltas(2,:),'o-')

n=48;
ts=linspace(8.25e-4,9.75e-4,n);
Deltas=zeros(2,n);
p.nfinal=250;
p.nradials=120;
p.lattice_symmetry='mm';
p.T=0.1*kB;
parfor i=1:n
    tic
    p1=p;
    p1.ts=ts(i);
    [Ds,~,~,~] = GKTH_self_consistency_iterate(p1,[D1,S2]);
    Deltas(:,i)=Ds;
    disp("ds iterate: "+i+" of "+n)
    toc
end
save("results\ss_bilayer\ds_Delta0_t_iter_zoom.mat",'ts','Deltas')
figure(3)
hold on
semilogx(ts,Deltas(1,:),'o-')
semilogx(ts,Deltas(2,:),'o-')

n=48;
ts=linspace(8.25e-4,9.75e-4,n);
Deltas=zeros(2,n);
p.nfinal=250;
p.nradials=120;
p.lattice_symmetry='mm';
p.T=0.1*kB;
parfor i=1:n
    tic
    p1=p;
    p1.ts=ts(i);
    [Ds,~,~,~] = GKTH_self_consistency_iterate(p1,[S1,D2]);
    Deltas(:,i)=Ds;
    disp("sd iterate: "+i+" of "+n)
    toc
end
save("results\ss_bilayer\sd_Delta0_t_iter_zoom.mat",'ts','Deltas')
figure(4)
hold on
semilogx(ts,Deltas(1,:),'o-')
semilogx(ts,Deltas(2,:),'o-')

%% Delta_0 calculation
n=48;
ts=logspace(-5,-1,n);
Deltas=zeros(2,n);
p.nfinal=250;
p.nradials=60;
p.lattice_symmetry='4mm';
p.T=0.1*kB;
parfor i=1:n
    tic
    p1=p;
    p1.ts=ts(i);
    [Ds,~,~] = GKTH_self_consistency_2S_taketurns(p1,[S1,S2]);
    Deltas(:,i)=Ds;
    disp("ss: "+i+" of "+n)
    toc
end
save("results\ss_bilayer\ss_Delta0_t.mat",'ts','Deltas')

n=48;
ts=logspace(-5,-1,n);
Deltas=zeros(2,n);
p.nfinal=250;
p.nradials=60;
p.lattice_symmetry='4mm';
p.T=0.1*kB;
parfor i=1:n
    tic
    p1=p;
    p1.ts=ts(i);
    [Ds,~,~] = GKTH_self_consistency_2S_taketurns(p1,[D1,D2]);
    Deltas(:,i)=Ds;
    disp("dd: "+i+" of "+n)
    toc
end
save("results\ss_bilayer\dd_Delta0_t.mat",'ts','Deltas')