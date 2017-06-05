% waveform_jweed.m %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script is to perform preliminary steps on looking at waveforms.sac
% after downloading them from jweed program
%
% what this script does:
% 1) Read SAC files
% 2) Apply Filters to seismograms
% 3) Choose to plot seismograms (filtered + unfiltered)
% 4) Pick or reject seismograms
% 5) Choose to record the selection of seismograms into an output file
%
% how to use this script:
% tweak the preliminary variables as needed
%
% Joshua Purba
% Canoe Reach Moment Tensor Project
% 6/1/2017
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear,clc
close all


%-------------------------Preliminary variables (change these as needed)
ievent=10;%Choose which EVENT to be used (ls from the directory of the event folder)
stnplot=100;%Choose number of STATIONS to analize and plot
iplot_single=0;%Choose to plot seismograms for one single station (unfiltered + filtered)
iplot_all=1;%Choose to plot seismograms for all stations 
iout=0; %Choose to write an output file or not
lowf=0.3;%Choose filter
highf=0.8;%Choose filter
xmin = 130; %p is p-arrival which around 150.0 sec
xmax = 250;

%-------Routine for organizing directory for events
eid_dir = dir('/home/joshuapurba/jweed_events/2*'); %find directory for events
eid_all = length([eid_dir.name]); %gives the length of char for all eids
eid_size = length([eid_dir(1).name]); %gives the length of char per eid
eid_num = eid_all/eid_size; %gives the number of events in such directory

%-------Routine for using SAC files
%for kk=1:eid_num; %loop for the entire events (kk is event #)
for kk=ievent; %loop for a specific event

%eid_sac is to get BHZ.sac for each event
eid_sac= dir(strcat('/home/joshuapurba/jweed_events/',eid_dir(kk).name,'/*BHZ*'));
cd(strcat('/home/joshuapurba/jweed_events/',eid_dir(kk).name)); %go to the directory of SAC
mm= length(eid_sac); %mm is the amount of BHZ.sac for each event

%-------Routine for sorting distance
for pp=1:mm %loop to get all distances from SAC files
    [sachdr,~] = load_sac(eid_sac(pp).name); %This is to read SAC heather and waveforms 
    stn_dist(pp)=(sachdr.dist); %save distance for every station
end
[stn_dist_sort,stn_idx] = sort(stn_dist'); %sort station in incremental epicentral distance

%-------Routine for stnplot correction (not that major)
if stnplot >= mm
    fprintf('stnplot is greater or equal to the number of stations which is %d\n',mm);
    fprintf('change stnplot into mm which is %d\n\n',mm);
    stnplot=mm;
else
    fprintf('stnplot is less than the number of stations which is %d\n\n',mm);
end

%-----loop for each SAC file (station) 
 for ll=1:stnplot
    [sachdr,data] = load_sac(eid_sac(stn_idx(ll)).name); %This is to read SAC heather and waveforms 

    %time used in seismogram, note that the +sachdr.delta means + 1 for index 0
    time = sachdr.b:sachdr.delta:sachdr.e + sachdr.delta; 
    if length(data) ~= length(time) %This is to check if the length of data=time
    time = sachdr.b:sachdr.delta:sachdr.e;
    end
    
    %finding p-arrival
    p = sachdr.t1; %p arrival is set to be 150 sec for the closest station
    p_ind = find(abs(time - p) == min(abs(time-p))); %finding P arrival indices
    p = p_ind*sachdr.delta; %The real P-arrival time
    
    %applying the pre-filter 
    offset = mean(data);
    data=data-offset; 

    %appying the filter
    %lowf=0.3;
    %highf=1; 
    dt=sachdr.delta; % need to check
    fs=1/dt;
    nyq=fs/2;
    [B,A]=butter(4,[lowf/nyq highf/nyq]);

    temp_filt = filtfilt(B,A,data);

    %----------------Routine for plotting one single station
    sachdr.knetwk=strtrim(sachdr.knetwk);
    sachdr.kstnm=strtrim(sachdr.kstnm);
    sachdr.kcmpnm=strtrim(sachdr.kcmpnm);
    if iplot_single 
    %ll = 1; mm =1; %##CHANGE THIS
   
    figure;
    %subplot(mm,2,ll+ll-1); %initially
    subplot(1,2,1);
    plot(time,data) %plot unfiltered seismogram
    line([p p],[-max(abs(data)) max(abs(data))],'Color','red','LineWidth',2); %plot p-arrival
    xlim([xmin xmax]);
    ylim([-max(abs(data)) max(abs(data))]);
    title(sprintf(' Unfiltered %s.%s.%s %0.2f km',sachdr.knetwk,sachdr.kstnm,sachdr.kcmpnm,sachdr.dist))

    %subplot(mm,2,ll+ll); %initially
    subplot(1,2,2);
    plot(time,temp_filt) %plot filtered seismogram
    line([p p],[-max(abs(temp_filt)) max(abs(temp_filt))],'Color','red','LineWidth',2); %plot p-arrival
    xlim([xmin xmax]);
    title(sprintf(' Filtered %s.%s.%s %0.2f km\nfilter: lowf=%.2fHz, highf=%.2fHz',sachdr.knetwk,sachdr.kstnm,sachdr.kcmpnm,sachdr.dist,lowf,highf))
    ylim([-max(abs(temp_filt)) max(abs(temp_filt))]);
    end
    
    %----------Routine for plotting all stations (filtered + unfiltered)
    if iplot_all  
    iplot_all_steps = (1+2*(ll-1)); %steps for plotting all seismograms
    
    figure(100);
    plot(time,iplot_all_steps+data/max(abs(data))) %plot unfiltered seismogram
    line([p p],[0 iplot_all_steps+1],'Color','red','LineWidth',2); %plot p-arrival
    title(sprintf('Unfiltered event year:%d day:%d hour:%d min:%d sec:%d \n lat:%.2f long:%.2f depth:%.2fkm\nchannel:%s',...
    sachdr.nzyear,sachdr.nzjday,sachdr.nzhour,sachdr.nzmin, sachdr.nzsec,sachdr.evla,sachdr.evlo,sachdr.evdp,sachdr.kcmpnm));
    xlim([xmin xmax]); xlabel(sprintf('Time windows in sec'));
    iplot_stn_text=sprintf('%s.%s %0.2f km',sachdr.knetwk,sachdr.kstnm,sachdr.dist);
    text(xmin,iplot_all_steps,iplot_stn_text,'HorizontalAlignment','right');%text for station
    set(gca,'YTickLabel',[]);
    hold on
    
    figure(101);
    plot(time,iplot_all_steps+temp_filt/max(abs(temp_filt))) %plot unfiltered seismogram
    line([p p],[0 iplot_all_steps+1],'Color','red','LineWidth',2); %plot p-arrival
    title(sprintf('Filtered event year:%d day:%d hour:%d min:%d sec:%d\nlat:%.2f long:%.2f depth:%.2fkm\nfilter: lowf=%.2fHz, highf=%.2fHz channel:%s',...
    sachdr.nzyear,sachdr.nzjday,sachdr.nzhour,sachdr.nzmin, sachdr.nzsec,sachdr.evla,sachdr.evlo,sachdr.evdp,lowf,highf,sachdr.kcmpnm));
    xlim([xmin xmax]); xlabel(sprintf('Time windows in sec'));
    text(xmin,iplot_all_steps,iplot_stn_text,'HorizontalAlignment','right')
    set(gca,'YTickLabel',[]);
    hold on
    end
 
    figure(102)
    plot(sachdr.evlo,sachdr.evla,'pentagram','color','r','linewidth',6)
    hold on
    plot(sachdr.stlo,sachdr.stla,'v','color','b','linewidth',6)
    hold on
    xlim([(sachdr.evlo)-5 (sachdr.evlo)+5]);
    ylim([(sachdr.evla)-5 (sachdr.evla)+5]);
    title(sprintf('Map of the event: year:%d day:%d hour:%d min:%d sec:%d \n lat:%.2f long:%.2f depth:%.2fkm',...
    sachdr.nzyear,sachdr.nzjday,sachdr.nzhour,sachdr.nzmin, sachdr.nzsec,sachdr.evla,sachdr.evlo,sachdr.evdp));
    xlabel(sprintf('longitude in degree'));
    ylabel(sprintf('latitude in degree'));
    iplot_stn_text2=sprintf('%s.%s\n',sachdr.knetwk,sachdr.kstnm);
    text(sachdr.stlo,sachdr.stla,iplot_stn_text2,'HorizontalAlignment','center')

      
    %---------Routine for keeping stations
    if iout
    keep=input(sprintf('No.%d keep this station? Y/N: ',ll),'s');
    if keep == 'Y' || keep =='y'
    ikeep(stn_idx(ll))=1; 
    ikeep_knetwk{stn_idx(ll)}=sachdr.knetwk;
    ikeep_kstnm{stn_idx(ll)}=sachdr.kstnm;
    ikeep_kcmpnm{stn_idx(ll)}=sachdr.kcmpnm;
    ikeep_dist(stn_idx(ll))=sachdr.dist;
    ikeep_stla(stn_idx(ll))=sachdr.stla;
    ikeep_stlo(stn_idx(ll))=sachdr.stlo;
    else
    ikeep(stn_idx(ll))=0;
    ikeep_knetwk{stn_idx(ll)}=[];
    ikeep_kstnm{stn_idx(ll)}=[];
    ikeep_kcmpnm{stn_idx(ll)}=[];
    ikeep_dist(stn_idx(ll))=NaN;
    ikeep_stla(stn_idx(ll))=NaN;
    ikeep_stlo(stn_idx(ll))=NaN;
    end
    end
  end
end

%-----------------------------------------Routine for writing .out file
if iout
fid=fopen(strcat(eid_dir(kk).name,'.out'),'w');

fprintf(fid,'output file for Canoe Reach Project \n\n');
fprintf(fid,'This is a summary for %s \nStation listed here are good and be used later:\n',eid_dir(kk).name');
fprintf(fid,'Event lat:%0.2f long:%0.2f \n\n\',sachdr.evla,sachdr.evlo);
fprintf('\n\n\nList of stations kept:\n\n')

ikeep_size = length(ikeep(find(ikeep~=0)));

for cc=1:mm
   if isempty(ikeep_knetwk{stn_idx(cc)})
   elseif ikeep(stn_idx(cc))==1; 
fprintf(fid,'%s.%s.%s %0.2fkm lat:%0.2f long:%0.2f\n',ikeep_knetwk{stn_idx(cc)},ikeep_kstnm{stn_idx(cc)},ikeep_kcmpnm{stn_idx(cc)},ikeep_dist(stn_idx(cc)),ikeep_stla(stn_idx(cc)),ikeep_stlo(stn_idx(cc)));
fprintf('%s.%s.%s %0.2fkm lat:%0.2f long:%0.2f\n',ikeep_knetwk{stn_idx(cc)},ikeep_kstnm{stn_idx(cc)},ikeep_kcmpnm{stn_idx(cc)},ikeep_dist(stn_idx(cc)),ikeep_stla(stn_idx(cc)),ikeep_stlo(stn_idx(cc)));
   
   end
end
fprintf(fid,'\nRatio (good)/(all) stations: %d/%d\n',ikeep_size,mm);
fprintf('\nRatio (good)/(all) stations: %d/%d\n',ikeep_size,mm);
fclose(fid);
end



%==========Below is to do test the amplitude to determine if such

%seismometer is good or not

%P_ind = find(abs(time - p) == min(abs(time-p))) %finding P arrival indices
%Pratio_aft = rms(data(P_ind+20:P_ind+500))
%Pration_bef =rms(data(P_ind-200:P_ind-20))
%plotting

% if iplot 
% figure(kk);
% subfigure(mm,1,ll);%mm=the number of stations,1=data,ll=seismometer no.ll
% plot(time,data)
% %xlim([147 169]);
% title('Filtered seismogram')
% 
% subfigure(mm,2,ll); %mm=the number of stations,2=filtered data,ll=seismometer no.ll
% plot(time,temp_filt)
% %xlim([147 169]);
% title('Filtered seismogram')
% end