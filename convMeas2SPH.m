%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% convMeas2SPH.m
% 
% Convert multichannel SOFA RIRs to Spherical harmonics (SH) SOFA  
%
% Author: Jung-Woo Choi (KAIST) and Franz Zotter (IEM, Graz)
% jwoo@kaist.ac.kr, zotter@iem.at

%% requirements
% (1) mtimesx library: for matrix multiplication over tensor data
%     https://www.mathworks.com/matlabcentral/fileexchange/
%     25977-mtimesx-fast-matrix-multiply-with-multi-dimensional-support
%
% (2) Spherical array processing toolbox by Archontis Politis 
%     https://github.com/polarch/Spherical-Array-Processing
%     https://github.com/polarch/Array-Response-Simulator
%     https://github.com/polarch/Spherical-Harmonic-Transform
%     (Archontis Politis, Microphone array processing for parametric 
%     spatial audio techniques, 2016,  
%     Doctoral Dissertation, Department of Signal Processing and Acoustics,
%     Aalto University, Finland)


clear 
close all

%% Parameters
%-  mode config
btrim = true; % trim data?
bplot = true; % plot IR?

ndat = 1024; % IR length to process & output 
tag = 'direct';
    if btrim, tag = ['_' tag '_trim' num2str(ndat)]; end
    
%- Data trimming config
Ttrim  = 0.7; % trimming size in seconds
nfadein   = 32;  % fade in length in samples
nfadeout  = 64;  % fade out length in samples
noffset = -16; % trim position offset in samples
    
%- path config
measname = 'measurement-no-obstacles'; %'measurement-wedges'; %'loc_test_full_run'; 

%- SPH encoding setting
N = 3; % encoding order
    nharm = (N+1)^2;  % num of SPH coefficients; 
sn = 20; % max eq gain (dB)
nenc = 1024; % length of encoding filter (samples)


%% define path & load data
datipath= ['.\SOFA_' measname '\'];
datopath= ['.\SOFA_SH_' measname '\'];
SOFAifn = [datipath measname]; 
SOFAofn = [datopath 'SH_' measname tag]; 
    if ~isfolder(datopath), mkdir(datopath); end

%- load sample data & shorten variable names 
Obj = SOFAload([SOFAifn '-a1.sofa']); % src, mic, dat

c0 = Obj.soundspeed;
rmic = Obj.micradius;
nmic = Obj.API.R;
nstep = Obj.nstep;
nspk = Obj.nspk;
nsrc = nspk*nstep;
narr = Obj.numarray;
fs   = Obj.Data.SamplingRate;

nkeep = round(Ttrim/c0*fs);  % win length to keep from direct

%% begin loop (array)
preVW = [0 0 0]; % listener view vector 
preUP = [0 0 0]; % listener up vector

for ia = 1:narr
    
%% load SOFA data
disp(['loading data:array' num2str(ia)]);
Obj = SOFAload([SOFAifn '-a' num2str(ia) '.sofa']); % src, mic, dat
disp('loading data complete');
ndat_org = size(Obj.Data.IR,3); % total num of data
if ~btrim, ndat = ndat_org; end

%- reshape info
arrpos = Obj.ListenerPosition(1,:); 
srcpos = Obj.SourcePosition;   % nsrc, 3

%% gen encoding mtx (Requires SHtoolbox)
%- update encoding filter if array orientation is different
if norm(preVW-Obj.ListenerView)~=0 || norm(preUP-Obj.ListenerUp)~=0
    disp('generating encoding filter')
    %- rotate the local mic coordinates to global
    micpos = rot_loc2glo(Obj.ReceiverPosition, Obj.ListenerView, Obj.ListenerUp);
    [~, encfilt] = arraySHTfiltersTheory_regLS_real(...
        rmic, micpos(:,1:2)/180*pi, N, nenc, fs, sn); %nm, mic, time
end

%% data trim
h = Obj.Data.IR;

if btrim 
    disp('trimming')
    nbegin = zeros(nsrc,1); % window beginning position

    %- position vector
    D = vecnorm(srcpos - arrpos,2,2); % src, 1

    %- gen truncation window (based on the geometric distance between SMA & src
    nbegin = floor(D/c0*fs)+noffset; % begin pass time
    trwin  = zeros(nsrc, ndat); % truncation window (src, time)
    for isrc =1:nsrc
        nbegin_ = max([1 nbegin(isrc)]); 
        if nbegin_ ~= 1
            nfadein_ = min([nbegin_-1 nfadein]);
            trwin(isrc, nbegin_+(-nfadein_:-1)) = ...
                fliplr(cos((1:nfadein_)/nfadein_*pi/2)).^2; % quater cosine win
        end
        nend_ = min([ndat nbegin_+nkeep-1]);
        if nend_ ~= ndat
            nfadeout_ = min([ndat-nend_ nfadeout]);
            trwin(isrc, nend_+(1:nfadeout_)) = ...
            cos((1:nfadeout_)/nfadeout_*pi/2).^2; % quater cosine win
        end
        trwin(isrc, nbegin_:nend_) = 1;
        trwin(isrc, (nend_+nfadeout_):end) = 0; 
    end

    %- multiply truncation window
    Obj.Data.IR = h(:,:,1:ndat).*permute(trwin,[1 3 2]); % src, mic, time
end

%% apply SPH encoder
disp('encoding')
%- take FFT of everything
Nfft = 2^(ceil(log2(size(encfilt,1) + ndat)));
Hw  = fft(Obj.Data.IR,Nfft,3);  % src, mic, freq 
EncFilt = fft(encfilt,Nfft,3);  % nm, nmic, freq (SH encoder)
predelay = size(encfilt,3)/2;   % delay by encoding filter
%- encoding 
Hnm = mtimesx(EncFilt, reshape(Hw, nsrc, nmic, Nfft),'t'); %  nm, src, freq
Hnm = reshape(Hnm,nharm,nsrc,Nfft); % nm, src, freq
hnm = ifft(Hnm,[],3);  % nm, src, time
hnm = circshift(hnm,[0 0 -predelay]);
hnm = permute(hnm(:,:,1:ndat),[2 1 3]); % src, nm, time (back to original convention)

%% Replace SOFA file and save
Listener = [ 'Encoded from Zylia ZM-1s into Spherical Harmonic Format using ACN and SN3D conventions. ' ...
               'SH Order:' num2str(N) ];
Obj.GLOBAL_ListenerDescription = Listener;

%- mic
Obj.ReceiverPosition = zeros(nharm,3);  
Obj.ReceiverPosition_Type = 'spherical harmonics';
Obj.ReceiverPosition_Units = 'ACN';   

%- IR
Obj.Data.IR = reshape(hnm, nsrc, nharm, []); % meas, nm, time
Obj.Data.Delay = zeros(nsrc,nharm,1);

%- Harmonics info
Obj = SOFAaddVariable(Obj,'Nenc','I',N); % encoding order

%- save
Obj = SOFAupdateDimensions(Obj);
disp(['Saving:  ' SOFAofn '-A' num2str(ia) ]);
Obj=SOFAsave([SOFAofn '-A' num2str(ia) '.sofa'], Obj);


%% plot data
if bplot
plotfun = @(x) 10*log10(abs(x));

figure(1)
set(gcf, 'position', [700 400 2150 930])
clf 
%- plot var
Hplot = squeeze(sum(h.^2,2))';
tax = (0:size(Hplot,1)-1)/fs*1000;
ylm = [0 40]; 
    
subplot(121)
imagesc(1:nsrc, tax, plotfun(Hplot)); % plot w ch, first step response (1:(Nspk*Nstep))/Nspk+1
if btrim
    hold on 
    winbound = (nbegin + [-nfadein 0 nkeep nfadeout])/fs*1000;
    hp = plot(1:nsrc, winbound);
    set(hp, {'color'}, num2cell(gray(4), 2));
end

title(['SMA ' num2str(ia) ' mic response'],'fontsize', 24); 
set([gca ylabel('time (ms)') xlabel('loudspeaker index')],'fontweight','bold','fontsize',20);
caxis([-60 0]+max(caxis))
ylim(ylm)
hb=colorbar; xlabel(hb,'(dB)')

subplot(122)
Hplot2 = squeeze(squeeze(hnm(:,1,:).^2))';
tax2 = (0:size(Hplot2,1)-1)/fs*1000;
    
imagesc(1:nsrc, tax2, plotfun(Hplot2)); % plot w ch, first step response (1:(Nspk*Nstep))/Nspk+1
title(['SMA ' num2str(ia) ' SPH response'],'fontsize',24); 
set([gca ylabel('time (ms)') xlabel('loudspeaker index')],'fontweight','bold','fontsize',20);
caxis([-60 0]+max(caxis))
ylim(ylm)
hb=colorbar; xlabel(hb,'(dB)')

drawnow

set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 11 9]);
print('-dpng','-r300',[SOFAofn '_IR_' num2str(ia) '_' num2str(N) '.png'])
saveas(gcf,[SOFAofn '_IR_' num2str(ia) '_' num2str(N) '.fig'])

%%
figure(2)
clf 
EDC = cumsum(squeeze(hnm(:,1,:).^2)','reverse'); % EDC of the zero-the order
EDCdB = 10*log10(EDC); 
EDCdBn = EDCdB-EDCdB(1,:);
tax3 = (0:size(EDCdBn,1)-1)/fs*1000;

imagesc(1:nsrc, tax3, EDCdBn); 
title(['Energy decay curve of array ' num2str(ia)]); ylabel('time (ms)'); xlabel('loudspeaker index');
caxis([-60 0]+max(caxis))

colorbar

end % end bplot

end % end ia


%% functions
function micpos = rot_loc2glo(micpos, vw, up)

vx = vw(:)'; vx = vx/norm(vx);
vz = up(:)'; vz = vz/norm(vz);
vy = cross(vz,vx); vy = vy/norm(vy); % gen local y axis

%- mic position in cartesian coordinates
[ex, ey, ez] = sph2cart(micpos(:,1)/180*pi,micpos(:,2)/180*pi,ones(size(micpos,1),1));

%- gen rotation matrix
R  = [vx; vy; vz];    % rotation matrix
gxyz = [ex ey ez]*R;  % transform mic positions to global coordinates

[maz, mel, ~] = cart2sph(gxyz(:,1), gxyz(:,2), gxyz(:,3));

micpos(:,1) = maz/pi*180; 
micpos(:,2) = mel/pi*180;
end