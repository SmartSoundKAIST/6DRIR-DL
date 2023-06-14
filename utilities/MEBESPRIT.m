clear 
close all

% Multi-array Eigenbeam ESPRIT code 
%
% For details of the algorithm, please see:
%
% J. -W. Choi, F. Zotter, B. Jo and J. -H. Yoo, 
% "Multiarray Eigenbeam-ESPRIT for 3D Sound Source Localization 
% With Multiple Spherical Microphone Arrays," 
% in IEEE/ACM Transactions on Audio, Speech, and Language Processing, 
% vol. 30, pp. 2310-2325, 2022, doi: 10.1109/TASLP.2022.3183929
%
% Author: Jung-Woo Choi (KAIST, South Korea) and Franz Zotter (IEM, Graz, Austria)


%% path config
addpath(genpath('.\lib')); % 

measname = 'measurement-no-obstacles'; %'measurement-wedges'; 

datipath = ['..\SOFA_SH_' measname '\']; % input data path (SH encoded SOFA)
datopath = ['..\MEB_ESPRIT_' measname '\']; %output data path
envpath  = ['..\SOFA_' measname '\'];   % environment file path
SOFAifn = [datipath 'SH_' measname '_direct_trim1024'];  % SOFA input file name 

%% setup
N= 3;        % max harmonics order to use
Nfft = 1024; % num of fft
ka = [1.8 2.2]; % frequency smoothing range
idxarr = [1:7]; % microphone arrays to use
nqa = 1;     % assumed num of sources

%- load environment file and shorten variable names
Obj = SOFAload([SOFAifn '-A1.sofa']); % src, mic, dat

c0 = Obj.soundspeed;
rmic = Obj.micradius;
nmic = Obj.API.R;
nstep = Obj.nstep;
nspk = Obj.nspk;
nsrc = nspk*nstep;
narr = Obj.numarray;
fs   = Obj.Data.SamplingRate;
Nenc = Obj.Nenc;

%- load room data
load([envpath measname '.mat']);
XJ = Env.jigpos; 
XR = Env.roomvertex;

%% Load data
disp('Loading data...');

%- build array structure
arr.n = numel(idxarr);
arr.xyz = zeros(numel(idxarr),3); % array pos
UV = zeros(numel(idxarr),3); VW = UV;

for ia = 1:numel(idxarr)
    
    Obj = SOFAload([SOFAifn '-a' num2str(idxarr(ia)) '.sofa']);
    Hnm{ia} = fft(Obj.Data.IR,Nfft,3);   %  meas, nharm, freq
    arr.xyz(ia,:) = Obj.ListenerPosition(1,:); % array position
    UV(ia,:) = Obj.ListenerUp;    % upvector (only for plotting)
    VW(ia,:) = Obj.ListenerView;  % view vector (only for plotting)
end
disp('Data loading complete');

%- source structure
src.xyz = Obj.SourcePosition; % src position
src.n = nsrc;

%% spectral config
fax = (0:Nfft-1)/Nfft*fs; fax = fax(1:(Nfft/2+1));
kax = 2*pi*fax/c0*rmic;

if numel(ka)>1
    kidx = find( ge(kax,ka(1)) & le(kax,ka(2))); 
else
    [~, kidx] = min(abs(kax-ka));
end

%% build cov mtx (with fsmoothing)
disp('Building covariances');

for iq = 1:src.n
    for ia1 = 1:arr.n
        for ia2 = 1:arr.n
            anm1 = squeeze(Hnm{ia1}(iq,:,kidx)); % nm, freq
            anm2 = squeeze(Hnm{ia2}(iq,:,kidx)); 
            Rt{ia1,ia2,iq} = (anm1*anm2')/numel(kidx); % nm, nm, Nfft, src
        end
    end
end
disp('Cov mtx generated')

%% MEB-ESPRIT
[C{1}, C{2}, C{3}] = multiplication_theorems_real(N);
M0 = [eye(N^2) zeros(N^2, 2*N+1)];

for iq = 1:nsrc % DoA estimation for each source
    iq
    R = Rt(:,:,iq);

    %- MEB-ESPRIT
    ESTm(iq) = mebESPRIT_randinit2022(R, M0, C, nqa, arr);
end

%% plot 
figure(1)
    clf
    clear pos

for ii=1:nstep*nspk
   prob_ = sum(ESTm(ii).prob.^2,2); % Qa,1
   prob(ii,:) = prob_;

   posall = ESTm(ii).xyz; 
   gtpos = src.xyz(ii,:); 
    
   [err(ii),idx(ii)] = min(sqrt(sum(abs(posall-gtpos).^2,2)));
   pos(ii,:) = posall(idx(ii),:);
end

erravg = rms(err); 
errstd = std(err);
[erravg errstd]

eplot = pos-src.xyz;
hqv= quiver3(src.xyz(:,1),src.xyz(:,2),src.xyz(:,3),eplot(:,1),eplot(:,2),eplot(:,3),0);
set(hqv,'linewidth',1.,'color','k');
title(['Position error (mean=' num2str(erravg*100, '%2.1f') ' cm, std.=' num2str(errstd*100, '%2.1f') ' cm)']);
xlabel('x (m)'); ylabel('y (m)');
axis equal
hold on 

plotconfig(arr.xyz, UV, VW,idxarr);
axis([-0.5    0.5   -1    1   -0.4    0.6]*1.6)

set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 11 9]*2);
print('-dpng','-r300','mebesprit_err.png')
saveas(gcf,'Measurement_config.fig')


function hnd = plotconfig(XA, UV, VW, idxarr)
%% plot config

%- origin & axis
scatter3(0,0,0,20,'mo','filled');
hold on
hnd.gax(1)=line([0 0.3],[0,0],[0,0],'color','b','linewidth',2);
hnd.gax(2)=line([0 0], [0,0.3],[0,0],'color','g','linewidth',2);
hnd.gax(3)=line([0 0], [0,0],[0,0.3],'color','r','linewidth',2);


%- SMA positions
[xp,yp,zp] = sphere(50);
sc = 0.054;
colors = [1 1 1 1 1 1 1 0.5 0.5];
for ii=1:size(XA,1)
    hnd.arr(ii) = surf(sc*xp+XA(ii,1),sc*yp+XA(ii,2),sc*zp+XA(ii,3), ones(size(xp))*colors(ii));
end
shading flat

% local axis
xloc = XA + VW*0.15;
zloc = XA + UV*0.15;
xcd = cat(3,XA,xloc); xcd = reshape(permute(xcd,[3,2,1]),6,[]);
zcd = cat(3,XA,zloc); zcd = reshape(permute(zcd,[3,2,1]),6,[]);
for ii=1:size(XA,1)
    hnd.arrax(1) = line(xcd(1:2,ii),xcd(3:4,ii),xcd(5:6,ii),'color','b','linewidth',2);
    hnd.arrax(2) = line(zcd(1:2,ii),zcd(3:4,ii),zcd(5:6,ii),'color','r','linewidth',2);

    if ii<8, zoff = -0.1; else zoff= 0.1; end
    hnd.arrtext(ii) = text(XA(ii,1), XA(ii,2), XA(ii,3)+zoff, [' A' num2str(idxarr(ii))]);
end
set(hnd.arrtext,'fontweight','bold','HorizontalAlignment','center')

axis equal
grid on
 
set([gca xlabel('x (m)') ylabel('y (m)')  zlabel('z (m)')],'fontweight','bold',...
 'fontsize',13)

axis equal
grid on
rotate3d on
camlight(135, 45)

camproj perspective
view(120,18)


end
