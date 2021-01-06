function outvec = aw_frequencySliding(data, freq_range, srate)
%
% outvec = aw_frequencySliding(data, freq_range, srate)
%
% measures instantaneous frequency by measuring the derivative of the phase
% angle timeseries according to methods of Cohen, M, J Neuro, 2014
% 
%       INPUTS: 
%               data =                 EEG data in EEGLAB format 
%                                           (channels x timepoints x trials)
%
%               freq_range =       vector of two integers (ex: [7 13]) 
%
%               srate =                 sampling rate of the data (ex: 500)
%                                               in Hz
%
%       OUTPUTS:
%               outvec =            matrix of calculated frequency sliding for the specified data
%                                           (channels x timepoints x
%                                           trials) ** NOTE: the timepoints
%                                           = input timepoints - 1
%
% written by annemarie wolff, awolf037@uottawa.ca
%% Determine the filter buffer for each edge of the filter
buffer_percent = .15;
buffer_hz = [freq_range(1)*buffer_percent, freq_range(2)*buffer_percent];

%% compute function
% initialize matrix
if ndims(data) == 3
    outvec = NaN([size(data,1), size(data,2)-1, size(data,3)]);
elseif ndims(data) == 2
    outvec = NaN([size(data,1), size(data,2)-1]);
else
    error('ERROR: data input dimensions are incorrect.');
end
           
tic;
for chani=1:size(data,1)
        % bandpass filter the data
        EEG_filt = bandpass(squeeze(data(chani,:,:)), [(freq_range(1)-buffer_hz(1)) (freq_range(2)+buffer_hz(2))], srate);
        
       %% take the hilbert tranform
        hilb_data = hilbert(EEG_filt);
        
       %%  compute frequency sliding
       if ndims(data) == 3
           % initialize matrix
           freqslide = NaN([size(data,2)-1 size(data,3)]);
            for triali=1:size(data,3)
                freqslide(:,triali) = srate*diff(unwrap(angle(hilb_data(:,triali))))/(2*pi);
            end
       else
                freqslide = srate*diff(unwrap(angle(hilb_data)))/(2*pi);
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
disp([ 'File completed frequency sliding in ' num2str(round(elapsed_time)) 's.']);
