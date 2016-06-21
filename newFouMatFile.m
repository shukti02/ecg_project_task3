function newFouMatFile(patient,idx)

matObj = matfile(patient);
filtMatObj = matfile('filteredLeads.mat');

if(nargin < 2)
    idx = [1 length(filtMatObj.I)];
end

featMatObj = matfile(strcat(patient, '_fou.mat'), 'Writable', true);
leads = {'I', 'II', 'III', 'aVF', 'aVL', 'aVR', 'V1', 'V2', 'V3', 'V4', 'V5', 'V6'};

beatpos = floor(1000* matObj.MSmodelsX);
beatPos_ms = beatpos(idx(1) <= beatpos & beatpos <= idx(end));
rrBord = [beatPos_ms(1:end-1); beatPos_ms(2:end)]';
rrBord(rrBord(:, 2) - rrBord(:, 1) > 3000, :) = [];
rrBord(rrBord(:, 2) - rrBord(:, 1) < 200, :) = [];
rrLen = double(rrBord(:, 2) - rrBord(:, 1) + 1);
rrCent = mean(rrBord, 2);
clear beatpos beatpos_ms
pt = floor(length(rrBord(:,2))/3);


for l = 1:length(leads)
    featMatObj.(strcat(leads{l}, '_fftMat'))(1,1:99) = 0;
    for i = 1 : 3
        signal = double(filtMatObj.(leads{l})(1,rrBord((i-1)*pt+1,1):rrBord((i)*pt+1,2)));
        rrBordi = rrBord((i-1)*pt+1:i*pt,:)-rrBord((i-1)*pt+1,1)+1;
        rrCenti = mean(rrBordi, 2);
        rrFft(signal,rrCenti,rrBordi,featMatObj,leads,l,i,pt);
    end
end

featMatObj.rrLoc = int32(idx(1) + rrCent' - 1);
featMatObj.rrLen = int32(rrLen');
saveArToFile(matObj,idx,featMatObj);
end


function rrFft(signal,rrCent,rrBord,featMatObj,leads,l,i,pt)
rrCell = cell(size(rrCent(:, 1)));
rrBord1 = rrBord(:,1); rrBord2 = rrBord(:,2);

    parfor rr = 1:size(rrCent, 1)
       rrIdx = rrBord1(rr):rrBord2(rr);
       rrCell{rr} = signal(rrIdx);  
    end
    clear signal
    fftCell = cellfun(@fft, rrCell, 'UniformOutput', false);
    clear rrCell

    % take frequencies 2 to 100 and store it into matrix fftMat
    fftMat = zeros(size(rrCent, 1), 99);
    for rr = 1:size(rrCent, 1)
        fftMat(rr, :) = abs(fftCell{rr}(2:100));
    end
    featMatObj.(strcat(leads{l}, '_fftMat'))((i-1)*pt+1:i*pt,:) = fftMat;
end

function saveArToFile(matObj,idx,featMatObj)
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
end