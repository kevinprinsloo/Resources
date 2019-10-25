%
%      Test GCFBv209 (for GCFBpack) 
%      Toshio IRINO
%      Created:   15 Sep 2005 (for v205)
%      Modified:   7 Apr 2006 (v206, Compensation of Group delay OutMidCrct)
%      Modified:  23 Dec 2006 (v207, 'dynamic' rather than 'time-varying')
%      Modified:  18 Mar 2007 (check. modifying title)
%      Modified:  25 Jul 2007 (GCresp)
%      Modified:   5 Aug 2007 (Elapsed time)
%      Modified:  26 Jun 2013
%      Modified:  25 Nov 2013 (v209)
%      Modified:  18 Apr 2015 (v210, include GCresp in GCFBv210_SetParam )
%      Modified:  13 May 2015 (v210, debug & comparison)
%
%
%clear

%%%% Stimuli : a simple pulse train %%%%
fs = 48000;
Tp = 10; % (ms) 100 Hz pulse train
Snd1 = [1, zeros(1,Tp*fs/1000-1)];
Snd = [];

for nnn = 1:10,
%for nnn = 1:3,
    Snd = [Snd1 Snd];
end;
Tsnd = length(Snd)/fs;
disp(['Duration of sound = ' num2str(Tsnd*1000) ' (ms)']);


SigSPLlist = [40:20:80]


ElapsedTime1 = [];
ErrorVersion_dc = [];
ErrorVersion_sc = [];

cnt = 0;
for nfig = 1:2
figure(nfig);

for nnn = 1:length(SigSPLlist)
SigSPL = SigSPLlist(nnn);
Snd =  Eqlz2MeddisHCLevel(Snd,SigSPL);

%%%% GCFB %%%%
GCparam = []; % reset all
GCparam.fs     = fs;
GCparam.NumCh  = 100;
GCparam.FRange = [100, 6000];

GCparam.OutMidCrct = 'No';
%GCparam.OutMidCrct = 'ELC';

if nfig == 1, GCparam.Ctrl = 'dynamic'; % used to be 'time-varying'
else    GCparam.Ctrl = 'static'; % or 'fixed' 
end;

tic
[cGCout, pGCout, GCparam210, GCresp210] = GCFBv210(Snd,GCparam);
tm210 = toc;

SwCompare = 1;
if SwCompare == 1
    cnt = cnt+1;
    tic
    [cGCout209, pGCout209, GCparam209, GCresp209] = GCFBv210(Snd,GCparam);
    tm209 = toc;
    ErrorVersion_dc(cnt) = sqrt(mean(mean((cGCout-cGCout209).^2)))/sqrt(mean(mean((cGCout209).^2)));
    ErrorVersion_sc(cnt) = sqrt(mean(mean((pGCout-pGCout209).^2)))/sqrt(mean(mean((pGCout209).^2)));
    ElapsedTime1(cnt,1:2) = [tm210, tm209];
end;

tm = toc;
disp(['Elapsed time is ' num2str(tm,4) ' (sec) = ' ...
      num2str(tm/Tsnd,4) ' times RealTime.']);

subplot(length(SigSPLlist),1,nnn)
imagesc(max(cGCout,0));
set(gca,'YDir','normal');
title(['GCFB control = "' GCparam.Ctrl '";  Signal Level = ' ...
       int2str(SigSPL) ' dB SPL']);

%GCparam210 
%GCresp210
drawnow
end;
end;
ErrorVersion_dc
ErrorVersion_sc
ElapsedTime1