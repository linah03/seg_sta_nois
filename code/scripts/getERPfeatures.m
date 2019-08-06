%This function loads EEGLAB *.set files and extracts features that are saved 
%to CSV file. The fnamefile input is a text file containing the filenames
%with paths to each of the *.set files to be loaded (one file per line).
%From a bash terminal these file name text files can be created as follows:
%for example...
%'find . -type f -name "*-SEG*.set" > derivatives/seg_face_noise/code/misc/fnames.txt
%The outfname input is the name of the file containing the output features.

function getERPfeatures(fnamefile,outfname)

fid=fopen(fnamefile,'r');
fnd=fread(fid);
fclose(fid);
cr_ind=find(fnd==10);
for i=1:length(cr_ind);
   if i==1;
      c_fname{i}=deblank(char(fnd(3:cr_ind(i))'));
      %EEG=eeg_checkset(EEG);
      %ALLEEG(i)=EEG;
   else
      c_fname{i}=deblank(char(fnd(cr_ind(i-1)+3:cr_ind(i))'));
      %EEG=eeg_checkset(EEG);
      %ALLEEG(i)=EEG;
   end
end


frontal_labs={'AF8','Fpz','AF7','F3','Fz','F4'};
central_labs={'C3','Cz','C4'};
posterior_labs={'P3','Pz','P4','PO7','PO8','Oz'};
temporal_labs={'FT8','TP8','FT7','TP7'};

fid=fopen(outfname,'w');
fprintf(fid,'%s',['fname,',...
                'frontal_area,central_area,posterior_area,temporal_area,', ...
                'frontal_peak,central_peak,posterior_peak,temporal_peak,', ...
                'frontal_peak_lat,central_peak_lat,posterior_peak_lat,temporal_peak_lat']);

for fi=1:length(c_fname);
    fprintf(fid,'\n%s,',c_fname{fi});
    EEG=pop_loadset('filename',c_fname{fi});
    EEG=eeg_checkset(EEG);
    
    for li=1:length(frontal_labs);
        frontal_inds(li)=find(strcmp({EEG.chanlocs.labels},frontal_labs{li}));
    end
    for li=1:length(central_labs);
        central_inds(li)=find(strcmp({EEG.chanlocs.labels},central_labs{li}));
    end
    for li=1:length(posterior_labs);
        posterior_inds(li)=find(strcmp({EEG.chanlocs.labels},posterior_labs{li}));
    end
    for li=1:length(temporal_labs);
        temporal_inds(li)=find(strcmp({EEG.chanlocs.labels},temporal_labs{li}));
    end
    
    ERP=mean(EEG.data,3);
    ERP_r(1,:)=mean(ERP(frontal_inds,:),1);
    ERP_r(2,:)=mean(ERP(central_inds,:),1);
    ERP_r(3,:)=mean(ERP(posterior_inds,:),1);
    ERP_r(4,:)=mean(ERP(temporal_inds,:),1);
    
    %0 to 800
    t_inds=find(EEG.times>0 & EEG.times<800);
    
    % measures
    % area [mean amplitude]
    area=mean(ERP_r(:,t_inds),2);
    
    % peakes
    %positive
    [peak_amp,peak_lat_ind]=max(ERP_r(:,t_inds),[],2);
    
    for ai=1:size(ERP_r,1);
        peak_lat(ai)=EEG.times(t_inds(peak_lat_ind(ai)));
    end
    
    
    % print current input file measures
    fprintf(fid,'%6.3f,',[area',peak_amp',peak_lat]);
    
end