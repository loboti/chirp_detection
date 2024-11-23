%% Program writes long-term EOD recordings to binary file and creates 
% log files in txt and in mat format with trial parameters

%% EDIT THE FOLLOWING 4 LINES BEFORE STARTING PROGRAM !!!
FishID = 12345678;  % fish ID
site = 'Tierhaus';
residentFishID = 'Q1,Q2,Q3,Q4,R1,R2,R3,R4, middle chamber';
neighborFishID = 'Q1,Q2,Q3,Q4,R1,R2,R3,R4, as plan';
Nside = 'Left'; %side chamber containing neighbor fish
Lside = 'Left'; %side chamber containing landmarks
tank = 'Q1,Q2,Q3,Q4,R1,R2,R3,R4'; %tanks associated to this recording
Temp = 24.6-25.7;             % water temperature
Conductivity = 190-210;    % water conductivity
pH = 8.5;              % water pH
DO = 7.0;              % dissolved oxygen (mg/l)

n_channels = 16;   % N of channels that get recorded, if you want to select specific channels on 2 DAQs, uncomment the relative lines (lines 25-42)
EODrec_amp = 'npi DPA-2FS';

%%
duration = 20; % max duration in seconds of data recorded into one bin file
imax = 168;     % max number of files to be written

d = daq.getDevices;
EODrec_DAQ = d.Model;   % write DAQ model to log file

s = daq.createSession('ni');
channels = 0:n_channels-1;  % data are written to channels AI1 and AI5
% s.addAnalogInputChannel(d.ID,channels, 'Voltage')
d1ch1=addAnalogInputChannel(s,'Dev1','ai0', 'Voltage')
d1ch2=addAnalogInputChannel(s,'Dev1','ai1', 'Voltage')
d1ch3=addAnalogInputChannel(s,'Dev1','ai2', 'Voltage')
d1ch4=addAnalogInputChannel(s,'Dev1','ai3', 'Voltage')
d1ch5=addAnalogInputChannel(s,'Dev1','ai4', 'Voltage')
d1ch6=addAnalogInputChannel(s,'Dev1','ai5', 'Voltage')
d1ch7=addAnalogInputChannel(s,'Dev1','ai6', 'Voltage')
d1ch8=addAnalogInputChannel(s,'Dev1','ai7', 'Voltage')

d1ch1.TerminalConfig = 'SingleEnded';
d1ch2.TerminalConfig = 'SingleEnded';
d1ch3.TerminalConfig = 'SingleEnded';
d1ch4.TerminalConfig = 'SingleEnded';
d1ch5.TerminalConfig = 'SingleEnded';
d1ch6.TerminalConfig = 'SingleEnded';
d1ch7.TerminalConfig = 'SingleEnded';
d1ch8.TerminalConfig = 'SingleEnded';

d2ch1=addAnalogInputChannel(s,'Dev2','ai0', 'Voltage')
d2ch2=addAnalogInputChannel(s,'Dev2','ai1', 'Voltage')
d2ch3=addAnalogInputChannel(s,'Dev2','ai2', 'Voltage')
d2ch4=addAnalogInputChannel(s,'Dev2','ai3', 'Voltage')
d2ch5=addAnalogInputChannel(s,'Dev2','ai4', 'Voltage')
d2ch6=addAnalogInputChannel(s,'Dev2','ai5', 'Voltage')
d2ch7=addAnalogInputChannel(s,'Dev2','ai6', 'Voltage')
d2ch8=addAnalogInputChannel(s,'Dev2','ai7', 'Voltage')

d2ch1.TerminalConfig = 'SingleEnded';
d2ch2.TerminalConfig = 'SingleEnded';
d2ch3.TerminalConfig = 'SingleEnded';
d2ch4.TerminalConfig = 'SingleEnded';
d2ch5.TerminalConfig = 'SingleEnded';
d2ch6.TerminalConfig = 'SingleEnded';
d2ch7.TerminalConfig = 'SingleEnded';
d2ch8.TerminalConfig = 'SingleEnded';



s.Rate = 20000;   % sample rate
Fs = s.Rate;

datapath = 'c:';
datapath = uigetdir(datapath,'Pick a directory to save acquired data');   % pick data directory 
datapath2 ='d:';
datapath2 = uigetdir(datapath2,'Pick a directory to be used in case the HD is full');   % pick data directory 


for i = 1:imax

            FileObj      = java.io.File(datapath);
            free_bytes   = FileObj.getFreeSpace
            total_bytes  = FileObj.getTotalSpace
            usable_bytes = FileObj.getUsableSpace


           if free_bytes <  1e+10
               datapath=datapath2;
           end
   
    cd(datapath);
    
    c = clock;
    logfilename = ['log_FishID' num2str(FishID) '_' num2str(c(1)) '_' num2str(c(2)) '_' num2str(c(3)) '_' num2str(c(4)) '_' num2str(c(5)) '.txt'];
    mat_logfilename = [logfilename(1:end-3) 'mat'];  % log variables in mat format
    datafilename = [logfilename(5:end-3) 'bin']    % data file name
    
    fid2 = fopen(logfilename, 'a');
    fprintf(fid2, '%s\r', [num2str(c(1)) '\' num2str(c(2)) '\' num2str(c(3)) ' ; ' num2str(c(4)) ':' num2str(c(5)) ':' num2str(c(6))]);
    fprintf(fid2, '%s%i\r', 'Resident Fish ID: ', residentFishID);
    fprintf(fid2, '%s%i\r', 'Neighbor Fish ID: ', neighborFishID);
    fprintf(fid2, '%s%i\r', 'Neighbor location: ', Nside);
    fprintf(fid2, '%s%i\r', 'Landmarks location: ', Lside);
    
    fprintf(fid2, '%s%s\r', 'Tank #: ', tank);
    fprintf(fid2, '%s%s\r', 'Data filename: ', datafilename);
    fprintf(fid2, '%s%i\r', 'Temperature: ', Temp);
    fprintf(fid2, '%s%i\r', 'Conductivity: ', Conductivity);
    fprintf(fid2, '%s%i\r', 'pH: ', pH);
    fprintf(fid2, '%s%i\r', 'DO: ', DO);
    fprintf(fid2, '%s%i\r', 'Sample rate: ', s.Rate);
    fprintf(fid2, '%s%i\r', '# channels: ', n_channels);
    fprintf(fid2, '%s%s\r', 'EOD recording amplifier: ', EODrec_amp);
    fprintf(fid2, '%s%s\r', 'EOD recording DAQ: ', EODrec_DAQ);
    fclose(fid2);
    

    save(mat_logfilename,'logfilename','datafilename','tank','residentFishID','neighborFishID','Nside','Lside','Temp',...
     'Conductivity','pH','DO','Fs','n_channels','EODrec_amp','EODrec_DAQ');
    fid1 = fopen(datafilename,'w');
    
   lh = s.addlistener('DataAvailable', @(src, event) logData(src, event, fid1));
   lh1 = s.addlistener('DataAvailable', @(src,event) plot(event.TimeStamps, event.Data));
    
    %% Acquire Data in the Background
    %
    % Acquire data continuously in a non-blocking mode.


    s.IsContinuous = true;
    s.startBackground;
   
   
    pause(duration);
    
    %% Stop the Continuous Acquisition and Close the Log File
    s.stop;
    delete(lh);
    delete(lh1);
    fclose(fid1);

end
%%
displayEndOfDemoMessage(datafilename)



function logData(src, evt, fid)
% Add the time stamp and the data values to data. To write data sequentially,
% transpose the matrix.
%   Copyright 2011 The MathWorks, Inc.

data = [evt.TimeStamps, evt.Data]' ;
fwrite(fid,data,'double');
end