%% loading data
done = 0;
load gradmat_IRT_spiral20shot.mat
input_waveform = gradmat(17,:);
ideal_waveform = gradmat(17,:);

%% add to path the required picoscope toolboxes
close all

% This section of code connects you to the TNMR .tnt file for read/write
ntnmr_app = actxserver('ntnmr.application');
ntnmr_data = actxserver('ntnmr.document');

% openfile MUST HAVE ABSOLUTE PATH
invoke(ntnmr_data,'OpenFile','D:\Jonathan\gradient_RL_lowfield\grad_measure_Y_noncartesian_NORF.tnt');

OFFSET = 4; % constant measurement offset between waveforms
DELAY = 0; % time delay between waveforms

%% loop: run TNMR, measure with picoscope (trap waveforms)

scan_duration= 3; % [s]
scan_dt = tnmr_timestring_to_float(invoke(ntnmr_data,'GetNMRParameter','Dwell Time')); % [s]
invoke(ntnmr_data,'SetNMRParameter','acqspiral', '6m');
acqspiral = tnmr_timestring_to_float(invoke(ntnmr_data,'GetNMRParameter','acqspiral')); % [s]

N_averages = 2;
average_pulse = zeros(3975,1);

pause(5);
for ii = 1:N_averages

    % set the gradient waveforms for both axes...put trig in front!:
    % this is in tnmr dt timescale
    trig = [0:0.0125:1,1:-0.0125:0];
    %ideal_p = [zeros(1,100),0:0.0125:1,ones(1,1000),1:-0.0125:0,zeros(1,100)]*80;
    ideal_p = input_waveform;
    test_waveform_length = length(ideal_p);
    ideal_p = [trig*50, zeros(1,1500), ideal_p];
    shapeP = mat2str(ideal_p);
    shapeP(1) = []; shapeP(end) = []; % trim '[' ']' characters
    invoke(ntnmr_data,'SetTable','shapeP', shapeP)
    invoke(ntnmr_data,'SetTable','shapeR', shapeP)
    % run the TNMR sequence
    disp(strcat(['Running TNMR seq (iter ',num2str(ii),')']))
    invoke(ntnmr_data,'ZG');

    % TODO: insert some update here (RL, etc.) of the gradient waveforms

    pause(scan_duration);
end
