%% loading data
done = 0;
load ideal_gradient_pulse.mat
load ideal_bipolar.mat
load gradmat_IRT_spiral20shot.mat
input_waveform = gradmat(17,:);
ideal_waveform = gradmat(17,:);

%% add to path the required picoscope toolboxes
close all
addpath(genpath("picoscope_instrument_control\"))
addpath(genpath("picoscope_support_toolbox\"))

% This section of code connects you to the TNMR .tnt file for read/write
ntnmr_app = actxserver('ntnmr.application');
ntnmr_data = actxserver('ntnmr.document');

% openfile MUST HAVE ABSOLUTE PATH
invoke(ntnmr_data,'OpenFile','D:\Jonathan\gradient_RL_lowfield\grad_measure_Y_noncartesian_NORF.tnt');

% CALIBRATION CONSTANTS - EMPIRICALLY DERIVED
TNMR_AMP_TO_PICO_MV = 12.7342; % conversion factor from TNMR grad waveform amplitude with
                               % waveform_amplitude = 100, TNMR A0= 14
                               % e.g. ~13 mV to waveform amplitude of 1 on scale
                               % [0 100]
RISE_EDGE_DELAY = 0.122*10^-3; %s the time in s when the trigger is picked up

OFFSET = 4; % constant measurement offset between waveforms
DELAY = 0; % time delay between waveforms

%% loop: run TNMR, measure with picoscope (trap waveforms)

scan_duration= 0.0; % [s]gb bbnn                           n
scan_dt = tnmr_timestring_to_float(invoke(ntnmr_data,'GetNMRParameter','Dwell Time')); % [s]
invoke(ntnmr_data,'SetNMRParameter','acqspiral', '6m');
acqspiral = tnmr_timestring_to_float(invoke(ntnmr_data,'GetNMRParameter','acqspiral')); % [s]

N_averages = 2;
average_pulse = zeros(3975,1);

% Initialize the Pico for measurement
PS4000A_Init_Block_Measure;
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

    % measure with the picoscope
    PS4000A_Block_Measure_NoOpenClose;

    pico_dt = timeIntervalNanoSeconds *10^-9; %s
    pulse_samples = (acqspiral+RISE_EDGE_DELAY)/pico_dt;

    pulse_out = chA(1:pulse_samples,:);
    average_pulse = average_pulse + pulse_out;

    % TODO: insert some update here (RL, etc.) of the gradient waveforms

    pause(scan_duration);
end
% Close the Picoscope
PS4000A_Close_Block_Measure;

average_pulse = average_pulse/N_averages;



%% PROCESSING THE DATA
close all
%TODO: Need to chek the details of this 
% Get the pulse to the correct time axis - registered to TNMR input pulse
t_pico = pico_dt:pico_dt:length(average_pulse)*pico_dt;
t_pico=t_pico+RISE_EDGE_DELAY; % actually measuring from trigger on

scan_dt = acqspiral/length(ideal_p);
t_tnmr = scan_dt:scan_dt:length(ideal_p)*scan_dt;
t_tnmr = t_tnmr+RISE_EDGE_DELAY;

average_pulse_t_tnmr = interp1(t_pico,average_pulse,t_tnmr);

first_trig_samp_t_tnmr = int16(RISE_EDGE_DELAY/scan_dt);

% figure, plot(t_tnmr, average_pulse_t_tnmr), 
ideal_p_from_trig = [ideal_p(first_trig_samp_t_tnmr+1:end),zeros(1,first_trig_samp_t_tnmr)];
% hold on, plot(t_tnmr,ideal_p_from_trig*TNMR_AMP_TO_PICO_MV)


measured_waveform = average_pulse_t_tnmr((length(average_pulse_t_tnmr)-test_waveform_length-first_trig_samp_t_tnmr):end-first_trig_samp_t_tnmr);
measured_waveform = measured_waveform(:,1:end-1)/TNMR_AMP_TO_PICO_MV;
measured_waveform = circshift(measured_waveform-OFFSET, -DELAY);
error = ideal_waveform - measured_waveform; 

%%
close all,
tiledlayout(1,2)
nexttile, plot(ideal_waveform)
hold on, plot(measured_waveform*1.15)

nexttile, plot(error), title('error (mV)')

savefig('current_design_output.fig')

save('current_measurement_data.mat', 'error', 'measured_waveform')
