clear,  daqreset; %close all, 
tic

cd('PATH_TO_YOUR_DATA')
%% select the channel to analyze (each tank has two channels)
n_channels=6;
activetank =1;      % indicate number of tank to analyze (this could be looped to extract data from all tanks), each tank correspond to 2 electrode pairs
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MANUAL=1; % select 1 for manual chirp detection, select 0 for automatic chirp detection
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs=20000;
freqrange = [0.5 1.5];  % frequency range to be displayed in spectrograms
nfft = 2^10; %2^14 gives the best frequency resolution but bad temporal one

th=10;

overlap = round(0.9*nfft); % um wieviel sollen Analysefenster ueberlappen? In Prozent
maxEOD1f = 1000;    % maximum frequency accepted as EODf1
minEOD1f = 700;     % minimum frequency accepted as EODf1
maxEOD2f = 1000;    % maximum frequency accepted as EODf2
minEOD2f = 600;     % minimum frequency accepted as EODf2
%% define segment to analyze

start_seg = 31 ; % starting point is set at this minute
end_seg   = 32 ; % end of analysis segment is set at this minute

%% open the file and search for the selected data range

[logfname, pathname] = uigetfile('*.mat','Pick a log file'); % verify that in the target folder the txt,mat and bin files correspond
logfullpathname = [pathname logfname];
load(logfullpathname);
datafullpathname = [pathname datafilename];
FID = fopen(datafullpathname,'r');

% fwrite(FID,[(start_seg*60*Fs):(end_seg*60*Fs-1)]);
% filename= ('C:\tmp\test.wav');
% audiowrite(filename,seg1,Fs);

%% load selected data segment and select channel to analyze: CHANNEL 1

currData = fread(FID,'double');
time = currData((1:(n_channels+1):end), 1);
EOD1  = currData(((2*activetank):(n_channels+1):end), 1);
trace1 = [time EOD1];
%crop1 = trace1((start_seg*60*Fs):(end_seg*60*Fs-1),:); % select data range of active channel to analyze
crop1 = trace1((start_seg*60*Fs-(Fs*60-1)):(start_seg*60*Fs),:);
seg1 = transpose(crop1(:, 2)); % should output column 2 values (EOD) within the time range defined by "crop"

%% load selected data segment and select channel to analyze: CHANNEL 2

EOD2  = currData(((2*activetank+1):(n_channels+1):end), 1);
trace2 = [time EOD2];
%crop2 = trace2((start_seg*60*Fs):(end_seg*60*Fs-1),:); % select data range of active channel to analyze
crop2 = trace2((start_seg*60*Fs-(Fs*60-1)):(start_seg*60*Fs),:);
seg2 = transpose(crop2(:, 2)); % should output column 2 values (EOD) within the time range defined by "crop"

%% calculate EOD freq
    
    EODsec1 = seg1(60*Fs-(Fs*60-1):60*Fs);
    [p1,f1] = pwelch(EODsec1,8192,4096,8192,Fs);      % extract EOD frequency, increae "nfft" and "window" in case of similar frequencies
    idx = find(p1==max(p1(find(f1<maxEOD1f & f1>minEOD1f))));
    EOD1f = f1(idx);
    
    EODsec2 = seg2(60*Fs-(Fs*60-1):60*Fs);
    [p2,f2] = pwelch(EODsec2, 8192,4096,8192, Fs);      % extract EOD frequency, increae "nfft" and "window" in case of similar frequencies
    idx = find(p2==max(p2(find(f2<maxEOD2f & f2>minEOD2f))));
    EOD2f = f2(idx);
    clear p1, clear f1, clear p2, clear f2, clear idx;
    
%% filter the two channels 
           band=5; % filter halfwidth: it determines the width of the filtered band and also the starting point of the chirp detection range
%           filt1=[EOD2f-band EOD2f+band];        %filter for EOD2f
%            filt2=[EOD1f-band EOD1f+band];        %filter for EOD1f    
%           filt3=[2*EOD2f-band 2*EOD2f+band];    %filter for ch2 second harmonic
%           filt4=[2*EOD1f-band 2*EOD1f+band];    %filter for ch1 second harmonic
%           seg1=bandstop(seg1,filt1,Fs); 
%            seg2=bandstop(seg2,filt2,Fs); 
 
filename= ('EOD.wav');
audiowrite(filename,[seg1+seg2],Fs);

    %% spectrogram plot - channel 1
figure();
ax1 = subplot(2,2,1);
title('spectrogram');

hold on
[spec1,f1,t1,p1]=spectrogram(seg1,blackmanharris(nfft),overlap,nfft,Fs,'yaxis','MinTHreshold',-80);
spectrogram(seg1,blackmanharris(nfft),overlap,nfft,Fs,'yaxis','MinTHreshold',-80);

ylim(freqrange);  % ordinate limits
set(gca, 'XTick', 0:10:60);          
set(gca, 'XTickLabel', 0:10:60);
hold off

hsp1 = get(gca, 'Position');  
set(gca,'Position', [hsp1(1), hsp1(2), hsp1(3), hsp1(4)]);

%% spectrogram plot - channel 2

ax3 = subplot(2,2,3);
title('spectrogram');
hold on

[spec2,f2,t2,p2]=spectrogram(seg2,blackmanharris(nfft),overlap,nfft,Fs,'yaxis','MinTHreshold',-80);
spectrogram(seg2,blackmanharris(nfft),overlap,nfft,Fs,'yaxis','MinTHreshold',-80);

ylim(freqrange);  % ordinate limits
set(gca, 'XTick', 0:10:60);          
set(gca, 'XTickLabel', 0:10:60);
hold off

hsp3 = get(gca, 'Position');  
set(gca, 'Position', [hsp1(1), hsp3(2), hsp1(3), hsp1(4)]);
    
%% add power plot for chirp detection using medfreq to extract the median freq from the spectrograms

ax2 =subplot(2,2,2);

fr1 = f1 > EOD1f+5 & f1< 2*EOD1f-100; %freq range in spec1
tr1 = t1 > 0 & t1 < 60;
m1 = medfreq(p1(fr1,tr1),f1(fr1));
% fr1a =  f1 > 2*EOD1f-15 & f1< 2*EOD1f+15
% h1a = medfreq(p1(fr1a,tr1),f1(fr1a),40);

fr2 = f2 > EOD2f+5 & f2< 2*EOD2f-100; %freq range in spec1
tr2 = t2 > 0 & t2 < 60;
m2 = medfreq(p2(fr2,tr2),f2(fr2));
% fr2a =  f2 > 2*EOD2f-15 & f2< 2*EOD2f+15
% h2a = medfreq(p2(fr2a,tr2),f2(fr2a),40);


%%
hold on

med1 = (m1-min(m1));
med2 = (m2-min(m2)); 

signal_1 = med1-med2;
signal_2 = med2-med1;

mw=(length(m1)/60)/1000*10; 

plot(t1(tr1),signal_1,'k');

[pks1,locs1,widths1]=findpeaks(signal_1,'MinPeakHeight',th,'MinPeakProminence',th, 'Annotate','extents','MinPeakDistance',2,'MinPeakWidth',mw/2); %,'MinPeakDistance',0.5,'MinPeakWidth',6
x_peaks1=t1(locs1);
chirps1=length(locs1); % net count of the chirps in channel 1 
 
ax21=scatter(x_peaks1,signal_1(locs1),'rv','MarkerFaceColor','r');
o1=isoutlier(signal_1(locs1),'mean');
ax31=scatter(x_peaks1(find(o1>0)),signal_1(find(o1>0)),'filled','MarkerFaceColor','b');
legend(ax21,['Chirp count =' num2str(chirps1)]);
hold off

title('resident chirp detection');
hsp2 = get(gca, 'Position');
set(gca, 'Position', [hsp2(1), hsp1(2), hsp1(3), hsp1(4)]);
ylabel = ('spec Power sum');

% add power plot for chirp detection
ax4 = subplot(2,2,4);

hold on
plot(t2(tr2),signal_2,'k');
% plot(t2(tr2),signal_ha_2,'b');

[pks2,locs2,widths2]=findpeaks(signal_2,'MinPeakHeight',th,'MinPeakProminence',th,'Annotate','extents','MinPeakDistance',2,'MinPeakWidth',mw/2); % 
x_peaks2=t2(locs2);
chirps2=length(locs2); % net count of the chirps in channel 2

ax41=scatter(x_peaks2,signal_2(locs2),'rv','MarkerFaceColor','r');
o1=isoutlier(signal_1(locs1),'mean');
ax31=scatter(x_peaks1(find(o1>0)),signal_1(find(o1>0)),'filled','MarkerFaceColor','b');

legend(ax41,['Chirp count =' num2str(chirps2)]);
hold off

title('neighbor chirp detection');
hsp4 = get(gca, 'Position');
set(gca, 'Position', [hsp2(1), hsp3(2), hsp1(3), hsp1(4)]);
ylabel = ('spec Power sum');

%% calculate EOD amplitude before each chirp
Alocs1=[];   %variable for signal amplitude 
Plocs1=[];   % variable for chirp max power


 
locs_sec1a=ceil(locs1/(length(m1)/60))*Fs - (Fs/1000)*3; % get sinewave fragments before chirps, using chirp timestamps. Removes 2 msec
locs_sec1b=locs_sec1a + (Fs/1000)*2; %create a second series of time stamps, delayed of 3 msec from the previous vector

for i=1:length(locs1);
 segs=seg1(locs_sec1a(i):(locs_sec1b(i)));   % extract signal segments


 
 [po,fo] = pwelch(segs); %determine freq and power of the analyzed segment

 
 Plocs1=[Plocs1, max(10*log10(po))]; % extract max power for analyzed segment
 Asegs=(max(segs)-min(segs))/2;      % calculate EOD amplitude from each segment
 Alocs1=[Alocs1, Asegs];             % append calculated values 
end

Alocs2=[];   %variable for signal amplitude
Plocs2=[];   % variable for chirp max power

 

 
locs_sec2a=ceil(locs2/(length(m2)/60))*Fs - (Fs/1000)*3; % get sinewave fragments before chirps, using chirp timestamps. Removes 2 msec
locs_sec2b=locs_sec2a + (Fs/1000)*2; %create a second series of time stamps, delayed of 3 msec from the previous vector

for i=1:length(locs2);
 segs2=seg2(locs_sec2a(i):(locs_sec2b(i)));   % extract signal segments


 
 [po2,fo2] = pwelch(segs2); %determine freq and power of the analyzed segment

 
 Plocs2=[Plocs2, max(10*log10(po2))]; %extract max power for analyzed segment
 Asegs2=(max(segs2)-min(segs2))/2; % calculate EOD amplitude from each segment
 Alocs2=[Alocs2, Asegs2];              %append calculated values 
end
%% convert variables
widths1=(widths1/195.9167)*1000; %convert amplitude in milliseconds
peaks1=m1(locs1)-EOD1f; % calculate real amplitudes of chirps taking them form the medfreq trace directly


widths2=(widths2/195.9167)*1000; %convert amplitude in milliseconds
peaks2=m2(locs2)-EOD2f; % calculate real amplitudes of chirps taking them form the medfreq trace directly


%% set all X axes equally scaled
allXLim = get([ax1, ax2, ax3, ax4], {'XLim'});
allXLim = cat(2, allXLim{:});
set([ax1, ax4], 'XLim', [min(allXLim), max(allXLim)]);
linkaxes([ax1,ax2,ax3,ax4],'x'); % this sets all X axes at the same scale

allYLim = get([ax1, ax3], {'YLim'});
allYLim = cat(2, allYLim{:});
set([ax1,ax3], 'YLim', [min(allYLim), max(allYLim)]);
linkaxes([ax1,ax3],'y'); % this sets all X axes at the same scale




toc