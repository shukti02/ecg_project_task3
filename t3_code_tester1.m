%tic;
clear all; clc;
patient = 'MG002';
matObj = matfile(patient);
filtMatObj = matfile('filteredLeads.mat');
% ecgFeat = matfile('ECGfeatures.mat');

%if(nargin < 2)
     idx = 1:length(filtMatObj.I);
%end

featMatObj = matfile(strcat(patient, '_fou.mat'), 'Writable', true);

%% save arrhythmia annotations to new file
 variables = fieldnames(matObj);
for v = 1:length(variables)
    if(regexp(variables{v}, 'arr\w*_x'));
        events_x = int32(1000*matObj.(variables{v}));
        events_len = int32(1000*matObj.(strcat(variables{v}(1:end-1), 'length')));
        events_sev = int32(1000*matObj.(strcat(variables{v}(1:end-1), 'severity')));
        
        eventIdx = idx(1) <= events_x & events_x <= idx(end); 
        
        featMatObj.(variables{v}) = events_x(eventIdx);
        featMatObj.(strcat(variables{v}(1:end-1), 'length')) = events_len(eventIdx);
        featMatObj.(strcat(variables{v}(1:end-1), 'severity')) = events_sev(eventIdx);
    end   
    
    if(regexp(variables{v}, 'GlucoseLevel'))
        gluc = matObj.GlucoseLevel;
        gluc_x = int32(1000*matObj.GlucoseLevelX);
        glucIdx = idx(1) <= gluc_x & gluc_x <= idx(end);
        
        featMatObj.GlucoseLevel = gluc(glucIdx);
        featMatObj.GlucoseLevelX = gluc_x(glucIdx);
    end
end

%% save fourier coefficients of rr intervals
%beatpos = ecgFeat.rPeaks;
beatpos = floor(1000* matObj.MSmodelsX);
beatPos_ms = beatpos(idx(1) <= beatpos & beatpos <= idx(end));
% beatPos_ms = int32(beatpos(...
%     idx(1) <= beatpos & beatpos <= idx(end)));

% Borders and Centers of RR intervals
rrBord = [beatPos_ms(1:end-1); beatPos_ms(2:end)]';
rrBord(rrBord(:, 2) - rrBord(:, 1) > 3000, :) = [];
rrBord(rrBord(:, 2) - rrBord(:, 1) < 200, :) = [];


rrLen = double(rrBord(:, 2) - rrBord(:, 1) + 1);
rrCent = mean(rrBord, 2);

rrBord1 = rrBord(:, 1);
rrBord2 = rrBord(:, 2);
pt = floor(length(rrBord(:,2))/3);
a = rrBord1(pt+1); b = rrBord1(2*pt+1);

leads = {'I', 'II', 'III', 'aVF', 'aVL', 'aVR', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6'};
% for l = 1:length(leads)
    tic
    % store rr intervals into cell array so to enable the use of cellfun
    rrCell = cell(size(rrCent(:, 1)));
%     signal = filtMatObj.(leads{1});
%     signal = double(signal);
    
    signal1 = double(filtMatObj.(leads{1})(1,1:rrBord2(pt+1))); signal2 = double(filtMatObj.(leads{1})(1,rrBord1(pt+1):rrBord2(2*pt+1)));
    signal3 = double(filtMatObj.(leads{1})(1,rrBord1(2*pt+1):end));
    
    parfor rr = 1:size(rrCent, 1)
        
        if (rr<=pt)
            rrIdx = rrBord1(rr):rrBord2(rr);
            rrCell{rr} = signal1(rrIdx);
        elseif (rr <= 2*pt)
            rrIdx = rrBord1(rr):rrBord2(rr);
            rrIdx = rrIdx-a+1;
            rrCell{rr} = signal2(rrIdx);
        else
            rrIdx = rrBord1(rr):rrBord2(rr);
            rrIdx = rrIdx-b+1;
            rrCell{rr} = signal3(rrIdx);
        end  
    end
  clear signal1 signal2 signal3
    T1 = toc
    tic
    fftCell = cellfun(@fft, rrCell, 'UniformOutput', false);
    clear rrCell

    % take frequencies 2 to 100 and store it into matrix fftMat
    fftMat = zeros(size(rrCent, 1), 99);
    for rr = 1:size(rrCent, 1)
        fftMat(rr, :) = abs(fftCell{rr}(2:100));
    end
    clear fftCell
    
    featMatObj.(strcat(leads{1}, '_fftMat')) = fftMat;
    T2 = toc
% end
featMatObj.rrLoc = int32(idx(1) + rrCent' - 1);
featMatObj.rrLen = int32(rrLen');
%toc;