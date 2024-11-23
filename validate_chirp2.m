%% manual detection for channel 1

clear, clc, close all; 
tt=1.2; % time threshold used to estimate chirp duration 
tt2=4;
figure1=1;
figure2=1;

ax4=[];

[logfname, pathname] = uigetfile('*.mat','Pick a log file'); %select a mat file from the folder containing the recording data
cd(pathname);
oldpath = path;
path('pathname',oldpath);
load(logfname);

locs2b=locs2;
locs20=[];
locs1b=locs1;
locs10=[];


if figure1==1;
h=figure(6); 

% scrsz = get(0,'ScreenSize');  % get screen size
% set(6,'Position',scrsz/2);

set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.6, 0.3, 0.4, 0.7]);
%set(gcf, 'Resize', 'off');
cla
rv=0;
cbreak=[];
restart=[];
undo=[];
r=[];
e=[];
zero=[]; 
ty=[];
 i=1;
duration1=[];
editA=[];

while i <= length(locs1) 

if restart==1, 
    locs1=locs1b;, i=1;, restart=[];, rv=0;,r=[];,e=[];,type1=zeros(1,length(locs1));locs10=[];,
end
if undo==1, 
    locs1(i-2)=locs1b(i-2);, undo=[];, rv=0;,r=[];,e=[]; locs10(i-2)=[];,i=i-2;
end
start=1;
    if cbreak==1, break; end
    
            if locs1(i)<=3 & length(locs1)>2% this prevents the first value to create plotting issues when too little
                locs1(i)=0;
                i=i+1;
                %uiresume(gcf)
            end
            
            if locs1(i) < 10 & locs1(i) > 2;
             start=locs1(i)-2;
            elseif locs1(i) > 10
             start=10;
            end
            
            if ~isempty(r)
               rv=get(r,'value');   % recall the previous value given to the checkbox
               if rv==1,
                locs1(i-1)=0;       % DELETES EVENT
               end
            end
            
            locs10=find(locs1==0);
            zero=zeros(1,length(locs10));
                        
            if locs1b(i) + 2*start < length(t1) & locs1b(i) - start > 0;
            timeseg = t1 > t1(locs1b(i) - start) & t1 < t1(locs1b(i) + 2*start);
            else
            timeseg = t1 > t1(locs1b(i) - 1) & t1 < t1(locs1b(i) + 1);
            end
             
            if locs1b(i) + 10*start < length(t1) &  locs1b(i) - start > 0;
            timeseg_spec = t1 > t1(locs1b(i) - start) & t1 < t1(locs1b(i) + 10*start);
            else
            timeseg_spec = t1 > t1(locs1b(i) - 1) & t1 < t1(locs1b(i) + 1);
            end
            
            medfreq1=m1(timeseg);
            medfreq2=m2(timeseg);
            l=1:50:500;
            lx=1;
            limit=start;
            if abs(medfreq1(1)-medfreq1(end)) > limit & abs((medfreq1(1) + medfreq1(end))/2 - EOD1f) > limit 
                  
                while abs(medfreq1(1)-medfreq1(end)) > limit  & abs((medfreq1(1)+medfreq1(end))/2 - EOD1f) > limit 
                     
                     if lx <= length(l)
                     limit0 = limit+l(lx);
                     lx=lx+1;
                     else limit0=100;
                     end
                     
                     if locs1(i) + limit0 < length(t1) & locs1(i) - limit0 > 0
                     timeseg = t1 > t1(locs1(i) - limit0) & t1 < t1(locs1(i) + limit0);
                     medfreq1=m1(timeseg);
                     limit=limit0;
                     elseif locs1(i) + limit0 > length(t1) & locs1(i) - limit0 > 0
                     timeseg = t1 > t1(locs1(i) - limit0) & t1 < t1(locs1(i) + (length(t1)- locs1(i)-1));
                     medfreq1=m1(timeseg);
                     elseif locs1(i) + limit0 < length(t1) & locs1(i) - limit0 < 0
                     timeseg = t1 > t1(locs1(i) -round(locs1(i)/2)) & t1 < t1(locs1(i) + limit0); %
                     medfreq1=m1(timeseg);   
                     break
                     end
                  end
            end
            
                            %% method2 (calculate duration)
                    m_rise=[m1(locs1(i)-1)];
                    t_rise=[t1(locs1(i)-1)];
                    adr=1;
                    add=1;
                    
                    while m_rise  > (min(medfreq1 + tt2)) %accumulate points until m_rise is close to baseline + threshold
                        if locs1(i)-1-adr > 0
                           m_rise = [m_rise, m1(locs1(i)-1-adr)];
                           t_rise = [t_rise, t1(locs1(i)-1-adr)];
                        elseif locs1(i)-adr > 0
                           m_rise = [m_rise, m1(locs1(i)-adr)];
                           t_rise = [t_rise, t1(locs1(i)-adr)];
                        else 
                           m_rise = [m_rise, m1(locs1(i))];
                           t_rise = [t_rise, t1(locs1(i))];

                           break
                        end
                        
                        adr = adr + 1;
                    end

                    m_decay=[m1(locs1(i))];
                    t_decay=[t1(locs1(i))];

                    while m_decay > (min(medfreq1 + tt2)) %accumulate points until m_decay is close to baseline + threshold

                        m_decay = [m_decay, m1(locs1(i)+add)];
                        t_decay = [t_decay, t1(locs1(i)+add)];
                        add = add + 1;
                    end

                    % discrete peak coordinates
                    m_peak=[fliplr(m_rise), m_decay];
                    t_peak=[fliplr(t_rise), t_decay];
%                        m_peak = unique(m_peak);
%                        t_peak = unique(t_peak);
                    time_peak=linspace(min(t_peak),max(t_peak),50); % time vector used to calculate the interpolated peak
                    peak=interp1(t_peak,m_peak,time_peak,'pchip'); % interpolated peak
                    Zt2=time_peak(find(peak > max([m_peak(1) m_peak(end)])));
                    
                    if ~isempty(time_peak) & ~isempty(Zt2)
                    duration3=(Zt2(end)-Zt2(1))*1000;
                    else
                    duration3=widths1(i);
                    end
                    %%
            
            ax1=subplot(4,3,[1 2 4 5 7 8],'Color',[0.4 0.4 0.4]);
            title('Chirp instant frequency (medfreq of EOD)')
            ax1.GridColor = [0.8, 0.8, 0.8];  % [R, G, B] 
            
            ax1.GridLineStyle = '--';
            hold on 
            plot(t1(timeseg),m1(timeseg),'w');
            plot(time_peak,peak,'w','LineWidth',2);
            plot(t2(timeseg),m2(timeseg),'y');
            scatter(x_peaks1(i),m1(locs1b(i)),50,'rv','MarkerFaceColor','r');
            xt=min(xlim(ax1))+0.005;
            yt=max(ylim(ax1))-10;
            
            hold off
            legend('\color{white} ch1','\color{yellow} ch2','\color{red} chirp','Location','northwest');
            box off
            grid on
            grid minor
            set(gca,'XMinorTick','on');
            xlabel('sec');
            txt1 = {'Measured Amplitude: ',num2str(peaks1(i)), 'Measured Duration:  ',num2str(widths1(i))};
            text(0.65,0.9,txt1, 'Units', 'Normalized','Color','w');
            text(0.65,0.8,'suggested values', 'Units', 'Normalized','Color','y');
            annotation(h,'textbox',[0.45 0.730 0.080 0.025],'Color',[1 1 0],'String','Duration','FitBoxToText','off','EdgeColor',[1 1 0]);
            annotation(h,'textbox',[0.45 0.760 0.080 0.025],'Color',[1 1 0],'String','Amplitude','FitBoxToText','off','EdgeColor',[1 1 0]);
            
             %% side plots

             
             ax5=subplot(4,3,6);
            
             fs1 = f1 > EOD1f-50 & f1< EOD1f+200; %freq range in sper
             hold on
             if sum(timeseg_spec)>1
             s1=surf(p1(fs1,timeseg_spec));
             s1.EdgeColor = 'none';
             view(2)
             end
            set(gca,'Visible','off')
            txt2 = {'ch1'};
            te2=text(93.5,12.4,txt2,'Color','w');
            scatter(10,0.5,'filled','r^');
            hold off
            
            ax4=subplot(4,3,9);
            
            fs2 = f2 > EOD2f-50 & f2< EOD2f+200; %freq range in sper
            if sum(timeseg_spec)>1
            s2=surf(p2(fs2,timeseg_spec));
            s2.EdgeColor = 'none';
            view(2)
            end
            txt3 = {'ch2'};
            te3=text(93.5,11.4,txt3,'Color','w');
            set(gca,'visible','off');
            
             %%
            
             ax3=subplot(4,3,3);
             title('Chirp power traces + err.')
             nfft=2^10;
             overlap=round(0.9*nfft);
             Fs=20000;
             [spec,f,t,p]=spectrogram(seg1,blackmanharris(nfft),overlap,nfft,Fs,'yaxis','MinTHreshold',-80);
             aspec=abs(spec);
             
             frange1=find(f1>(EOD1f-10) & f1<(EOD1f+10));
             frange2=find(f1>(EOD1f+20) & f1<(EOD1f+500));
             
             basetrace = sum(aspec(frange1,:),1);
             basetrace_seg = basetrace(timeseg);
             chirptrace = sum(aspec(frange2,:),1);
             chirptrace_seg = chirptrace(timeseg);
             rejection=[];
             if length(chirptrace_seg)>1 & length(basetrace_seg) > 1
%              hold on
%              plot(basetrace_seg)
%              plot(chirptrace_seg)
%              peak=10; 
%              scatter(peak,0,'rv','MarkerFaceColor','r');
%              dif=chirptrace_seg-basetrace_seg;
%              dif=dif-min(dif);
%              plot(dif,'--')
%              legend('baseline', 'signal', 'chirp','diff') ;
%              hold off
              
                t_seg=t(timeseg); % chirp time vector
                ti_seg=linspace(t_seg(1),t_seg(end),100); %interpolate chirp time vector
                chirptrace_i_seg=interp1(t_seg,chirptrace_seg,ti_seg,'pchip'); %interpolated chirp vector
                basetrace_i_seg=interp1(t_seg,basetrace_seg,ti_seg,'pchip');  %interpolated baseline vector

                idxs=((find(chirptrace_i_seg > basetrace_i_seg))); % peak indexes, these can be used to compare basetrace and chirptrace in terms of averages
                idxz=((find(chirptrace_i_seg < basetrace_i_seg))); % peak indexes
                
                
                if ~isempty(idxs)
                    rejection=[];
                    diff_i = chirptrace_i_seg- basetrace_i_seg;

                    hold on
                    scatter(ti_seg(idxs),diff_i(idxs),10,'r','filled');
                    scatter(ti_seg(idxz),diff_i(idxz),10,'b','filled');
                    scatter(x_peaks1(i),0,'rv','MarkerFaceColor','r');
                    legend('chirp', 'base EOD', 'peak','Location','northeastoutside'); 
                    box off
                    
                    % 1st condition: chirp timestamp included in the chirp spectrogram time window
                    if x_peaks1(i) < ti_seg(idxs(1)) & x_peaks1(i) > ti_seg(idxs(end)) 
                    rejection=[rejection,1];
                    end
                    % 2nd condition: mean peak is higher than mean gap
                    if mean(chirptrace_i_seg(idxs)) < mean(chirptrace_i_seg(idxz))         
                    rejection=[rejection,2];
                    end

                    % 3rd condition: baseline excursion (max-min) 
                    if abs(max(diff_i(idxz))-min(diff_i(idxz))) > abs(max(diff_i(idxs))-min(diff_i(idxs)))
                        rejection=[rejection,3];
                    end
                    % 4th condition: comparison of absolute means
                    if abs(mean(chirptrace_i_seg(idxs))) < abs(mean(chirptrace_i_seg(idxz)))
                        rejection=[rejection,4];
                    end     
                    
                    bm=islocalmin(basetrace_i_seg(idxs));
                    cm=islocalmax(chirptrace_i_seg(idxs));
                    if sum(imdilate(bm,ones(3))).*(imdilate(cm,ones(3))) < 1
                    rejection=[rejection,5];
                    end 
                    
                    % 6th condition: peaks and gaps not neighbors
                    if abs(find(basetrace_i_seg==(min(basetrace_i_seg(20:end-40))))) - abs(find(chirptrace_i_seg==(max(chirptrace_i_seg(20:end-40))))) > 5
                    rejection=[rejection,6];
                    end 
                    % 7th condition: average chirp power 
                    if mean(chirptrace_i_seg(idxs)) < 5
                        rejection=[rejection,7];
                    end
                    % 8th condition: baseline length
                    if length(diff_i(idxz))<40
                        rejection=[rejection,8];
                    end
                else rejection=[rejection,9];
                end
             else rejection=[rejection,8];
             end
                if sum(rejection) >0
                  txt3 = {num2str(rejection)};
                  text(0.05,0.95,txt3, 'Units', 'Normalized','Color','r','LineWidth',2); 
                end
              
              
              
              
              
             %%
            ax2=subplot(4,3,[10 12],'Color',[0.4 0.4 0.4]); 
            ax2.GridColor = [0.8, 0.8, 0.8];  % [R, G, B] 
            grid on
            locs10=find(locs1==0);
            zero=zeros(1,length(locs10));
            
            hold on
            plot(t1,signal_1,'y');
            s3=scatter(t1(locs1b(i)),pks1(i),50,'filled','rv');
            s4=scatter(t1(locs1b(locs10)),zero,'filled','b');
            
            hold off
            legend('\color{yellow} ch1','\color{red} chrip','\color{blue} deleted');
            txt = ['Chirp # ' num2str(i) '/' num2str(length(locs1))];
            txto = ['total chirps = ' num2str(length(locs1b)) ' + ' num2str(length(locs2b))];
            text(0.05,0.9,txt,'units','normalized','Color','w');
            text(0.05,0.8,txto,'units','normalized','Color','w');
            ylim([0 150]);
            if t1((locs1(i))) < 20.0000;
                xlim([0 20.0000]);
            elseif t1((locs1(i))) < 40.0000 & t1((locs1(i))) > 20.0000;
                xlim([20.0000 40.0000]);
            elseif t1((locs1(i))) < 60.0000 & t1((locs1(i))) > 40.0000;
                xlim([40.0000 60.0000]);
            end
             xlabel('sec');
            %%
            warning('off','all');
            
            edit1=0; edit2=0;
            widths1(i)=duration3;
            ed=uicontrol('style','Edit','string',widths1(i),'Units','normalized','Position',[0.535 0.730 0.080 0.025],'callback','uicontrol(ed), widths1(i-1)=str2num(ed.String); edit1=1;');
            ea=uicontrol('style','Edit','string',max(m1(timeseg))-EOD1f,'Units','normalized','Position',[0.535 0.760 0.080 0.025],'callback','uicontrol(ea), peaks1(i-1)=str2num(ea.String); edit2=1;');
           
            if edit1==0

                if ~isempty(duration3);
                ed.String=duration3;
                widths1(i)=duration3;
                trace_m=m1(timeseg);
                amplitude1=max(m_peak)-EOD1f;  %mean(m2(timeseg)
                ea.String=amplitude1;
                end
            end
            
            if edit2==0
                peaks1(i)=max(m_peak)-EOD1f; %AMPLITUDE CALCULATION
                
            end
            
            
            
             
             
            if ~isempty(r)
               rv=get(r,'value');   % recall the previous value given to the checkbox
               if rv==1,
                locs1(i-1)=0; % DELETES EVENT
               locs10=find(locs1==0);
               end
            end
            a = uicontrol('Style','pushbutton','position',[10,10,60,50],'FontWeight','bold','String','Accept','Callback','uiresume(gcf)');
            e = uicontrol('Style','pushbutton','position',[940,10,60,25],'FontWeight','bold','String','END','Callback','cbreak=1; uiresume(gcf); close');  
            reset = uicontrol('Style','pushbutton','position',[940,130,60,25],'String','reset','Callback','restart=1; uiresume(gcf);');  
            transfer = uicontrol('Style','pushbutton','position',[940,160,60,25],'String','transfer','Callback','Plocs2=[Plocs2,Plocs1(i-1)]; Alocs2=[Alocs2,Alocs2(i-1)]; locs2=[locs2,locs1(i-1)]; sort(locs2);locs2b=locs2; x_peaks2=[x_peaks2,x_peaks1(i-1)]; sort(x_peaks2);pks2=[pks2,pks1(i-1)]; sort(pks2); peaks2=[peaks2,peaks1(i-1)]; sort(peaks2) ;widths2=[widths2,widths1(i-1)]; sort(widths2); locs1(i-1)=0; locs10=find(locs1==0)  ; uiresume(gcf);');
            undo = uicontrol('Style','pushbutton','position',[940,190,60,25],'String','undo','Callback','undo=1; uiresume(gcf);');  
            %% 
            
            k = uicontrol('Style','pushbutton','value',rv,'position',[10,60,60,50],'ForegroundColor',[1 0 0],'FontWeight','bold','String','Reject','Callback','locs1(i-1)=0;locs10=find(locs1==0);uiresume(gcf);');
            
            
            %sk = uicontrol('Style','pushbutton','position',[940,40,60,25],'String','skip+10','Callback',' i=i+9;  type1(i-10:i-1)= val;  uiresume(gcf);');
            df = uicontrol('Style','pushbutton','position',[940,40,60,25],'String','del+5','Callback',' locs1(i-1:i+3)= 0; i=i+4; uiresume(gcf);'); 
            dk = uicontrol('Style','pushbutton','position',[940,70,60,25],'String','del+10','Callback',' locs1(i-1:i+8)= 0; i=i+9;  uiresume(gcf);'); 
            dall = uicontrol('Style','pushbutton','position',[940,100,60,25],'String','del all','Callback',' locs10=find(locs1); cbreak=1; uiresume(gcf)'); 
            %%
            i=i+1;
            uiwait(gcf);
            
            if ishandle(ax1) & ishandle(ax2); % if the figure is still open, refresh the plots
            cla(ax1);
            cla(ax2);
            cla(ax4);
            cla(ax5);
            cla(ax3);
            end
            

          
end

%%
Alocs1(locs10)=[];

Plocs1(locs10)=[];
locs1(locs10)=[];
locs1(locs1==0)=[];
x_peaks1(locs10)=[];

peaks1(locs10)=[];
pks1(locs10)=[];
chirps1=length(locs1);

widths1(locs10)=[];
ICI1=diff(x_peaks1);
chirps1=length(locs1); % net count of the chirps in channel 1 
freq1=chirps1/60;
%%

disp(sprintf('channel 1 processed, now processing channel 2 '));
%pause
close ;
end

if figure2==1;
%%
h=figure(6);
% scrsz = get(0,'ScreenSize');  % get screen size
% set(6,'Position',scrsz);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.6, 0.4, 0.4, 0.6]);
%set(gcf, 'Resize', 'off');
cla
rv=0;
cbreak=[];
restart=[];
r=[];
e=[];
zero=[]; 
duration2=[];
 i=1;


while i <= length(locs2)

if restart==1, 
    locs2=locs2b;, i=1;, restart=[];, rv=0;,r=[];,e=[];,type2=zeros(1,length(locs2));locs20=[];,
end
if undo==1, 
    locs2(i-2)=locs2b(i-2);, undo=[];, rv=0;,r=[];,e=[]; locs20(i-2)=[];,i=i-2;
end

    if cbreak==1, break; end
    
            if locs2(i)<=3 & length(locs2)>2 % this prevents the first value to create plotting issues when too little
                locs2(i)=0;
                i=i+1;
            end
            
            if locs2(i) < 10;
             start=locs2(i)-1;
            else
             start=10;
            end
            
            if ~isempty(r)
               rv=get(r,'value');   % recall the previous value given to the checkbox
               if rv==1,
                locs2(i-1)=0;       % DELETES EVENT
               end
            end


            
            locs20=find(locs2==0);
            zero=zeros(1,length(locs20));
            
            if locs2b(i) + 2*start < length(t2) & locs2b(i) - start > 0;
            timeseg = t2 > t2(locs2b(i) - start) & t2 < t2(locs2b(i) + 2*start);
            else 
            timeseg = t2 > t2(locs2b(i) - 1) & t2 < t2(locs2b(i) + 1);
            end    
                
                
            if locs2b(i) + 10*start < length(t2) & locs2b(i) - start > 0;
            timeseg_spec = t2 > t2(locs2b(i) - start) & t2 < t2(locs2b(i) + 10*start);
            else
            timeseg_spec = t2 > t2(locs2b(i) - 1) & t2 < t2(locs2b(i) + 1);
            end
            
            if sum(timeseg_spec)<=1;
                
                uiresume(gcf)
            end
           %
            
           
            medfreq1=m1(timeseg);
            medfreq2=m2(timeseg);

            l=1:50:500;
            lx=1;
            limit=start;
            if abs(medfreq2(1)-medfreq2(end)) > limit & abs((medfreq2(1)+medfreq2(end))/2 - EOD2f) > limit 
                  while abs(medfreq2(1)-medfreq2(end)) > limit  & abs((medfreq2(1)+medfreq2(end))/2 - EOD2f) > limit % 
                     if lx <= length(l)
                     limit0 = limit+l(lx);
                     lx=lx+1;
                     else limit0=100;
                     end
                     if locs2(i) + limit0 < length(t2) & locs2(i) - limit0 > 0
                     timeseg = t2 > t2(locs2(i) - limit0) & t2 < t2(locs2(i) + limit0);
                     medfreq2=m2(timeseg);
                     limit=limit0;
                     elseif locs2(i) + limit0 > length(t2) & locs2(i) - limit0 > 0
                     timeseg = t2 > t2(locs2(i) -limit0) & t2 < t2(locs2(i) + (length(t2)- locs2(i)-1));
                     medfreq2=m2(timeseg);
                     elseif locs2(i) + limit0 < length(t2) & locs2(i) - limit0 < 0
                     timeseg = t2 > t2(locs2(i) -round(locs2(i)/2)) & t2 < t2(locs2(i) + limit0);
                     medfreq2=m2(timeseg);
                     break
                     end
                  end
            end
            % method2 (calculate duration)
            
                    m_rise=[m2(locs2(i)-1)];
                    t_rise=[t2(locs2(i)-1)];
                    adr=1;
                    add=1;
                    


                    
                     while m_rise  > (min(medfreq2 + tt2)) %accumulate points until m_rise is close to baseline + threshold
                        if locs2(i)-1-adr > 0
                           m_rise = [m_rise, m2(locs2(i)-1-adr)];
                           t_rise = [t_rise, t2(locs2(i)-1-adr)];
                        elseif locs2(i)-adr > 0
                           m_rise = [m_rise, m2(locs2(i)-adr)];
                           t_rise = [t_rise, t2(locs2(i)-adr)];
                        else 
                           m_rise = [m_rise, m2(locs2(i))];
                           t_rise = [t_rise, t2(locs2(i))];

                           break
                        end
                        
                        adr = adr + 1;
                    end
                    
                    
                    
                    
                    
                    
                    m_decay=[m2(locs2(i))];
                    t_decay=[t2(locs2(i))];
                    

                    while m_decay > (min(medfreq2 + tt2)) %accumulate points until m_decay is close to baseline + threshold
                        if locs2(i)+add < length(m2)
                        m_decay = [m_decay, m2(locs2(i)+add)];
                        t_decay = [t_decay, t2(locs2(i)+add)];
                        else 
                        m_decay = [m_decay, m2(locs2(i))];
                        t_decay = [t_decay, t2(locs2(i))];
                        break
                        end
                        add = add + 1;
                    end


                    % discrete peak coordinates
                    m_peak=[fliplr(m_rise), m_decay];
                    t_peak=[fliplr(t_rise), t_decay];
                       m_peak = unique(m_peak,'stable');
                       t_peak = unique(t_peak,'stable');
                       
                       if length(t_peak)>length(m_peak)
                           t_peak(length(m_peak)+1:end)=[];
                       end
                       if length(m_peak)>length(t_peak)
                           m_peak(length(t_peak)+1:end)=[];
                       end 
                       
                    time_peak=linspace(min(t_peak),max(t_peak),50); % time vector used to calculate the interpolated peak
                    if length(t_peak)==1 & length(m_peak)==1
                        t_peak(2)=t_peak(1)+1;
                        m_peak(2)=m_peak(1)+1;
                    end
                    
                    peak=interp1(t_peak,m_peak,time_peak,'pchip'); % interpolated peak
                    Zt2=time_peak(find(peak > max([m_peak(1) m_peak(end)])));
                    
                    if ~isempty(time_peak) & ~isempty(Zt2);
                    duration3=(Zt2(end)-Zt2(1))*1000;
                    else duration3=widths2(i);
                    end
           
         
           
            ax1=subplot(4,3,[1 2 4 5 7 8],'Color',[0.4 0.4 0.4]);
            title('Chirp instant frequency (medfreq of EOD)')
            ax1.GridColor = [0.8, 0.8, 0.8];  % [R, G, B] 
            
            ax1.GridLineStyle = '--';
            hold on 
            plot(t1(timeseg),m1(timeseg),'w');
            plot(t2(timeseg),m2(timeseg),'y');
            plot(time_peak,peak,'y','LineWidth',2);
            scatter(x_peaks2(i),m2(locs2b(i)),50,'rv','MarkerFaceColor','r');

            %timeseg = t2 > t2(locs2b(i) - start) & t2 < t2(locs2b(i) + 2*start);
            xt=min(xlim(ax1))+0.005;
            yt=max(ylim(ax1))-10;
            
            hold off
            legend('\color{white} ch1','\color{yellow} ch2','\color{red} chirp','Location','northwest');
            box off
            
            grid on
            grid minor
            set(gca,'XMinorTick','on');
            xlabel('sec');
            txt1 = {'Measured Amplitude: ',num2str(peaks2(i)), 'Measured Duration:  ',num2str(widths2(i))};
            text(0.65,0.9,txt1, 'Units', 'Normalized','Color','w');
            text(0.65,0.8,'suggested values', 'Units', 'Normalized','Color','y');
            annotation(h,'textbox',[0.45 0.730 0.080 0.025],'Color',[1 1 0],'String','Duration','FitBoxToText','off','EdgeColor',[1 1 0]);
            annotation(h,'textbox',[0.45 0.760 0.080 0.025],'Color',[1 1 0],'String','Amplitude','FitBoxToText','off','EdgeColor',[1 1 0]);
            % side plots

             ax5=subplot(4,3,6);
             %axes('Position',[.7 .6 .2 .2])
        if sum(timeseg_spec) >1;
             fs1 = f1 > EOD1f-50 & f1< EOD1f+200; %freq range in sper
            
             s1=surf(p1(fs1,timeseg_spec));
             s1.EdgeColor = 'none';
             view(2)

            set(gca,'Visible','off')
            txt2 = {'ch1'};
            te2=text(93.5,12.4,txt2,'Color','w');
            
            
            ax4=subplot(4,3,9);
            
            fs2 = f2 > EOD2f-50 & f2< EOD2f+200; %freq range in sper
            hold on
            s2=surf(p2(fs2,timeseg_spec));
            s2.EdgeColor = 'none';
            
            view(2)

            txt3 = {'ch2'};
            te3=text(93.5,11.4,txt3,'Color','w');
            scatter(10,0.5,'filled','r^');
            hold off
            set(gca,'visible','off');
            
        else uiresume(gcf);
        end  
            
%             ax3=subplot(4,3,3,'Color',[0.4 0.4 0.4]);

             ax3=subplot(4,3,3);
             title('Chirp power traces + err.')
             nfft=2^10;
             overlap=round(0.9*nfft);
             Fs=20000;
             [spec,f,t,p]=spectrogram(seg2,blackmanharris(nfft),overlap,nfft,Fs,'yaxis','MinTHreshold',-80);
             aspec=abs(spec);
             
             frange1=find(f2>(EOD2f-10) & f2<(EOD2f+10));
             frange2=find(f2>(EOD2f+20) & f2<(EOD2f+500));
             
             basetrace = sum(aspec(frange1,:),1);
             basetrace_seg = basetrace(timeseg);
             chirptrace = sum(aspec(frange2,:),1);
             chirptrace_seg = chirptrace(timeseg);
             rejection=[];
             
             if length(chirptrace_seg)>1 & length(basetrace_seg) > 1
%              hold on
%              plot(basetrace_seg)
%              plot(chirptrace_seg)
%              
%              peak=10; 
%              scatter(peak,0,'rv','MarkerFaceColor','r');
%              dif=chirptrace_seg-basetrace_seg;
%              plot(dif,'--')
%              legend('baseline', 'signal', 'chirp','dif') ;
%              hold off
                
                t_seg=t(timeseg); % chirp time vector
                ti_seg=linspace(t_seg(1),t_seg(end),100); %interpolate chirp time vector
                chirptrace_i_seg=interp1(t_seg,chirptrace_seg,ti_seg,'pchip'); %interpolated chirp vector
                basetrace_i_seg=interp1(t_seg,basetrace_seg,ti_seg,'pchip');  %interpolated baseline vector

                idxs=((find(chirptrace_i_seg > basetrace_i_seg))); % peak indexes, these can be used to compare basetrace and chirptrace in terms of averages
                idxz=((find(chirptrace_i_seg < basetrace_i_seg))); % peak indexes
                
                rejection=[];
                if ~isempty(idxs)

                    diff_i = chirptrace_i_seg- basetrace_i_seg;

                    hold on
                    scatter(ti_seg(idxs),diff_i(idxs),10,'r','filled');
                    scatter(ti_seg(idxz),diff_i(idxz),10,'b','filled');
                    scatter(x_peaks2(i),0,'rv','MarkerFaceColor','r');
                    legend('chirp', 'base EOD', 'peak','Location','northeastoutside');
                    box off
                    
                    % 1st condition: chirp timestamp included in the chirp spectrogram time window
                    
                    if x_peaks2(i) < ti_seg(idxs(1)) & x_peaks2(i) > ti_seg(idxs(end)) 
                    rejection=[rejection,1];
                    
                    end
                    % 2nd condition: mean peak is higher than mean gap
                    if mean(chirptrace_i_seg(idxs)) - mean(chirptrace_i_seg(idxz)) < 4        
                    rejection=[rejection,2];
                    end

                    % 3rd condition: baseline excursion (max-min) 
                    if abs(max(diff_i(idxz))-min(diff_i(idxz))) > abs(max(diff_i(idxs))-min(diff_i(idxs)))
                        rejection=[rejection,3];
                    end
                    % 4th condition: comparison of absolute means
                    if abs(mean(chirptrace_i_seg(idxs))) < abs(mean(chirptrace_i_seg(idxz)))
                        rejection=[rejection,4];
                    end     
                    
                    %5th condition: peaks and gaps not aligned
                    bm=islocalmin(basetrace_i_seg(idxs));
                    cm=islocalmax(chirptrace_i_seg(idxs));
                    if sum(imdilate(bm,ones(3))).*(imdilate(cm,ones(3))) < 1
                    rejection=[rejection,5];
                    end 
                    
                    % 6th condition: peaks and gaps not neighbors
                   if abs(find(basetrace_i_seg==(min(basetrace_i_seg(20:end-40))))) - abs(find(chirptrace_i_seg==(max(chirptrace_i_seg(20:end-40))))) > 5
                   rejection=[rejection,6];
                   end
                   
                   % 7th condition: average chirp power
                   if mean(chirptrace_i_seg(idxs)) < 5
                       rejection=[rejection,7];
                   end
                   % 8th condition: baseline length
                   if length(diff_i(idxz))<40
                       rejection=[rejection,8];
                   end
                else rejection=[rejection,9];
                end
             else rejection=[rejection,8];
             end
                    if sum(rejection) > 0
                      txt3 = {num2str(rejection)};
                      text(0.05,0.95,txt3, 'Units', 'Normalized','Color','r','LineWidth',2); 
                    end
             
            %
         
%

            ax2=subplot(4,3,[10 12],'Color',[0.4 0.4 0.4]); 
            ax2.GridColor = [0.8, 0.8, 0.8];  % [R, G, B] 
            grid on
            locs20=find(locs2==0);
            zero=zeros(1,length(locs20));
            
            hold on
            plot(t2,signal_2,'y');
            s3=scatter(t2(locs2b(i)),pks2(i),50,'filled','rv');
            s4=scatter(t2(locs2b(locs20)),zero,'filled','b');
            s5=scatter(t1(locs1),signal_1(locs1),20,'filled','w');            
            hold off
            legend('\color{yellow} ch2','\color{red} chrip','\color{blue} deleted','\color{white} ch1');
            txt = ['Chirp # ' num2str(i) '/' num2str(length(locs2))];
            txto = ['total chirps = ' num2str(length(locs1b)) ' + ' num2str(length(locs2b))];
            text(0.05,0.9,txt,'units','normalized','Color','w');
            text(0.05,0.8,txto,'units','normalized','Color','w');
            ylim([0 150]);
            if t2((locs2(i))) < 20.0000;
                xlim([0 20.0000]);
            elseif t2((locs2(i))) < 40.0000 & t2((locs2(i))) > 20.0000;
                xlim([20.0000 40.0000]);
            elseif t2((locs2(i))) < 60.0000 & t2((locs2(i))) > 40.0000;
                xlim([40.0000 60.0000]);
            end

             xlabel('sec');
         
            warning('off','all');
            
            edit1=0;edit2=0;
            ed=uicontrol('style','Edit','string',widths2(i),'Units','normalized','Position',[0.535 0.730 0.080 0.025],'callback','uicontrol(ed), widths2(i-1)=str2num(ed.String); edit1=1;');
            ea=uicontrol('style','Edit','string',max(m2(timeseg))-EOD2f,'Units','normalized','Position',[0.535 0.760 0.080 0.025],'callback','uicontrol(ea), peaks2(i-1)=str2num(ea.String);edit2=1;');
            

            if edit1==0
               
                
               if ~isempty(duration3);
                ed.String=duration3;
                widths2(i)=duration3;
                trace_m=m1(timeseg);
                amplitude1=max(m_peak)-EOD2f;
                ea.String=amplitude1;
                end
            end
            
            if edit2==0
                peaks2(i)=max(m_peak)-EOD2f;
            end
            
            %r = uicontrol('Style','checkbox','value',rv,'position',[90,10,50,20],'String','Reject','Callback','rv=1;');
             
             
            if ~isempty(r)
               rv=get(r,'value');   % recall the previous value given to the checkbox
               if rv==1,
                locs2(i-1)=0; % DELETES EVENT
               locs20=find(locs2==0);
               end
            end
            a = uicontrol('Style','pushbutton','position',[10,10,60,50],'FontWeight','bold','String','Accept','Callback','uiresume(gcf)');
            e = uicontrol('Style','pushbutton','position',[940,10,60,25],'FontWeight','bold','String','END','Callback','cbreak=1; uiresume(gcf); close');  %on=max(length(locs2)); uiresume(gcf)
            reset = uicontrol('Style','pushbutton','position',[940,130,60,25],'String','reset','Callback','restart=1; uiresume(gcf);');  %on=max(length(locs2)); uiresume(gcf)
            transfer = uicontrol('Style','pushbutton','position',[940,160,60,25],'String','transfer','Callback','Alocs1=[Alocs1,Alocs2(i-1)];Plocs1=[Plocs1,Plocs2(i-1)]; locs1=[locs1,locs2(i-1)]; locs2(i-1)=0; chirps1=length(locs1); freq1=chirps1/60; x_peaks1=[x_peaks1,x_peaks2(i-1)]; sort(x_peaks1);  ICI1=diff(x_peaks1); pks1=[pks1,pks2(i-1)]; peaks1=[peaks1,peaks2(i-1)]; widths1=[widths1,widths2(i-1)]; locs20=find(locs2==0) ; uiresume(gcf);');
            undo = uicontrol('Style','pushbutton','position',[940,190,60,25],'String','undo','Callback','undo=1; uiresume(gcf);');  
            k = uicontrol('Style','pushbutton','value',rv,'position',[10,60,60,50],'ForegroundColor',[1 0 0],'FontWeight','bold','String','Reject','Callback','locs2(i-1)=0;locs20=find(locs2==0);uiresume(gcf);');


%            sk = uicontrol('Style','pushbutton','position',[940,40,60,25],'String','skip+10','Callback',' i=i+9;  type2(i-10:i-1)= val;  uiresume(gcf);');
            df = uicontrol('Style','pushbutton','position',[940,40,60,25],'String','del+5','Callback',' locs2(i-1:i+3)= 0; i=i+4;  uiresume(gcf);'); 

            dk = uicontrol('Style','pushbutton','position',[940,70,60,25],'String','del+10','Callback',' locs2(i-1:i+8)= 0; i=i+9;  uiresume(gcf);'); 
            dall = uicontrol('Style','pushbutton','position',[940,100,60,25],'String','del all','Callback',' locs20=find(locs2); cbreak=1; uiresume(gcf)'); 
          
            i=i+1;
            uiwait(gcf);
            
            if ishandle(ax1) & ishandle(ax2) & ishandle(ax4) & ishandle(ax5); % if the figure is still open, refresh the plots
            cla(ax1);
            cla(ax2);
            cla(ax4);
            cla(ax5);
            cla(ax3);
            end
          
end


%%


Alocs2(locs20)=[];

Plocs2(locs20)=[];
locs2(locs20)=[];
x_peaks2(locs20)=[];
widths2(locs20)=[];
peaks2(locs20)=[];
pks2(locs20)=[];
MSlocs2(locs20)=[];
ICI2=diff(x_peaks2);
chirps2=length(locs2); % net count of the chirps in channel 1 
freq2=chirps2/60;


close ;
end



%%
clear seg1;
clear seg2;
clear spec;
clear aspec;
%%
validated_log=[logfname(1:end-4),'_checked']
save(validated_log);

