function outvec = aw_frequencyPower(data, freq_range, srate)
%
% outvec = aw_frequencyPower(data, freq_range, srate)
%
% measures instantaneous power by measuring square of the magnitude of the 
% analytic signal according to methods of Cohen, M, J Neuro, 2014
% 
%       INPUTS: 
%               data =                 EEG data in EEGLAB format 
%                                           (channels x timepoints x trials)
%
%               freq_range =       vector of two integers (ex: [7 13])
%                                           buffer of 15% for lower and
%                                           upper bounds will be added for
%                                           edge effects
%
%               srate =                 sampling rate of the data (ex: 500)
%                                               in Hz
%
%       OUTPUTS:
%               outvec =            matrix of calculated instantaneous
%                                           power for the specified data with the mean taken over all the trials
%                                           (channels x timepoints x trials) 
%
% written by annemarie wolff, awolf037@uottawa.ca
%% Determine the filter buffer for each edge of the filter
buffer_percent = .15;
buffer_hz = [freq_range(1)*buffer_percent, freq_range(2)*buffer_percent];

%% compute function
% initialize matrix
outvec = NaN(size(data));

tic;
for chani=1:size(data,1)
        % bandpass filter the data
        EEG_filt = bandpass(squeeze(data(chani,:,:)), [(freq_range(1)-buffer_hz(1)) (freq_range(2)+buffer_hz(2))], srate);
        
       %% take the hilbert transform
        hilb_data = hilbert(EEG_filt);
        
      %%  compute power sliding
       if ndims(data) == 3
           % initialize matrix
           freqslide = NaN([size(data,2)-1 size(data,3)]);
            for triali=1:size(data,3)
                freqslide(:,triali) = (abs(hilb_data(:,triali))).^2;
            end
       else
                freqslide = (abs(hilb_data)).^2;
       end
        clear hilb_data  
        
       %% apply a median filter, then take the mean over all the trials
       if ndims(data) == 3       
    %         outvec(chani,:) = mean(medfilt1(freqslide, 10, [], 2,'omitnan'),2);
            outvec(chani,:,:) = medfilt1(freqslide, 10, [], 2,'omitnan');
       else
             outvec(chani,:) = medfilt1(freqslide, 10, [], 2,'omitnan');
       end
        clear freqslide
          
end

elapsed_time = toc;
disp([ 'File completed frequency power in ' num2str(round(elapsed_time,1)) 's.']);
