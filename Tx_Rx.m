%% reset all
rng ('default'); clear ; close all ;

pkg load signal;
pkg load communications;

%% Part_1_Tx_UPNRZ

% Generating_random_stream_of_bits

tx_bit_stream_1 = randi([0, 1], 1,10000); % 10,000 random bits


% Line coding the stream of bits
V_levels_1 = [0 1.2]; % Voltage levels for '0'and '1'
line_coded_stream_1 = V_levels_1(tx_bit_stream_1 + 1); % Assign voltage levels based on bit values
Counter_1=1;
Counter_2=100;
tx_Uni_polar_Non_Return_Zero_signal=zeros(1,1000000);
for Counter_3= 1:length(line_coded_stream_1)
    for Counter_4 =Counter_1:Counter_2
        tx_Uni_polar_Non_Return_Zero_signal(Counter_4)=line_coded_stream_1(Counter_3);
    end
    Counter_1=Counter_1+100;
    Counter_2=Counter_2+100;
end

Counter_11=1;
Counter_22=100;
tx_bit_stream_11=zeros(1,1000000);
for Counter_33= 1:length(tx_bit_stream_1)
    for Counter_44 =Counter_11:Counter_22
        tx_bit_stream_11(Counter_44)=tx_bit_stream_1(Counter_33);
    end
    Counter_11=Counter_11+100;
    Counter_22=Counter_22+100;
end


% Time_Vector
Bits_duration =0.01;
t_1 = 0:Bits_duration:(length(tx_Uni_polar_Non_Return_Zero_signal) * Bits_duration) - Bits_duration;


% Time-Domain Plot
figure (11)
subplot(2,1,1)
plot(t_1,tx_Uni_polar_Non_Return_Zero_signal);
xlabel('Time');
ylabel('Voltage');
title('Uni-polar Non Return to Zero (NRZ) Line Coding Example');
xlim([0 5]);
ylim([-0.5 1.5]);
grid on;

subplot(2,1,2)
plot(t_1,tx_bit_stream_11);
xlabel('Time');
ylabel('Value');
title('The Code');
xlim([0 5]);
ylim([-0.5 1.5]);
grid on;

%Frequency Domain Plot
Fs=1/Bits_duration;
T=10000;
df=1/T;
N=ceil(T/Bits_duration);
if (rem(N,2)==0)
    freq=(-0.5*Fs):df:(0.5*Fs)-df; %even function
else
    freq=-(0.5*Fs-0.5*df):df:(0.5*Fs-0.5*df); %odd function
end
freq_spectrum_1 = abs(fftshift(fft(tx_Uni_polar_Non_Return_Zero_signal))).^2/N;
figure (9)
subplot(2,3,1)
plot(freq,freq_spectrum_1);
xlabel('Frequency');
ylabel('Magnitude squared');
title('Uni-polar Non Return to Zero');
grid on;



% Eye_diagram_plot_using_builtin_function
eyediagram(tx_Uni_polar_Non_Return_Zero_signal,300,3);
title('Eye-Diagram for Uni-polar Non Return to Zero');


% Transmit the signal
rx_Uni_polar_Non_Return_Zero_signal= tx_Uni_polar_Non_Return_Zero_signal;

% Sweep on the value of sigma and calculate BER for each value
sigma_values_1 = linspace(0, max(V_levels_1), 10);
ber_values_1 = zeros(1, length(sigma_values_1));
rx_bit_stream_1 = zeros (1,length(tx_bit_stream_1));
for i = 1:length(sigma_values_1)
    % Add noise to the received signal
    noise_1 = sigma_values_1(i) * randn(1, length(rx_Uni_polar_Non_Return_Zero_signal));
    rx_Uni_polar_Non_Return_Zero_signal_with_noise = rx_Uni_polar_Non_Return_Zero_signal+ noise_1;

    for index = 50:100:length(rx_Uni_polar_Non_Return_Zero_signal_with_noise)

        rx_bit_stream_1((index+50)/100) = (rx_Uni_polar_Non_Return_Zero_signal_with_noise(index )> 0.6);

    end

    % Compare received bit stream with transmitted bit stream and count errors
    num_errors_1 = sum(rx_bit_stream_1 ~= tx_bit_stream_1);
    % Calculate BER for current value of sigma
    ber_values_1(i) = num_errors_1 / length(tx_bit_stream_1);
end

% Plot BER vs. sigma
figure (8)
semilogy(sigma_values_1, ber_values_1);
xlabel('Noise RMS Value (sigma)');
ylabel('Bit Error Rate (BER)');
title('BER vs. Noise RMS Value');
grid on;


%% Part_1_PNRZ

% Generating_random_stream_of_bits

tx_bit_stream_2 = randi([0, 1], 1, 10000); % 10,000 random bits

V_levels_2 = [-1.2 1.2]; % Voltage levels for '0'and '1'
line_coded_stream_2 = V_levels_2(tx_bit_stream_2 + 1); % Assign voltage levels based on bit values
Counter_1=1;
Counter_2=100;
tx_Polar_Non_Return_Zero_signal=zeros(1,1000000);
for Counter_3= 1:length(line_coded_stream_2)
    for Counter_4=Counter_1:Counter_2
        tx_Polar_Non_Return_Zero_signal(Counter_4)=line_coded_stream_2(Counter_3);
    end
    Counter_1=Counter_1+100;
    Counter_2=Counter_2+100;
end

Counter_11=1;
Counter_22=100;
tx_bit_stream_22=zeros(1,1000000);
for Counter_33= 1:length(tx_bit_stream_2)
    for Counter_44 =Counter_11:Counter_22
        tx_bit_stream_22(Counter_44)=tx_bit_stream_2(Counter_33);
    end
    Counter_11=Counter_11+100;
    Counter_22=Counter_22+100;
end

% Time_Vector
Bits_duration =0.01;
t_2 = 0:Bits_duration:(length(tx_Polar_Non_Return_Zero_signal) * Bits_duration) - Bits_duration;


% Time-Domain Plot
figure (12)
subplot(2,1,1)
plot(t_2,tx_Polar_Non_Return_Zero_signal);
xlabel('Time');
ylabel('Voltage');
title('polar Non Return to Zero (NRZ) Line Coding Example');
ylim([-1.5 1.5]);
xlim([0 5]);
grid on;

subplot(2,1,2)
plot(t_2,tx_bit_stream_22);
xlabel('Time');
ylabel('Value');
title('The Code');
xlim([0 5]);
ylim([-0.5 1.5]);
grid on;

%Frequency Domain Plot
Fs=1/Bits_duration;
T=10000;
df=1/T;
N=ceil(T/Bits_duration);
if (rem(N,2)==0)
    freq=(-0.5*Fs):df:(0.5*Fs)-df; %even function
else
    freq=-(0.5*Fs-0.5*df):df:(0.5*Fs-0.5*df); %odd function
end
freq_spectrum_2 = abs(fftshift(fft(tx_Polar_Non_Return_Zero_signal))).^2/N;
figure (9)
subplot(2,3,2)
plot(freq,freq_spectrum_2);
xlabel('Frequency');
ylabel('Magnitude squared');
title('Polar Non Return to Zero');
grid on;


% Eye_diagram_plot_using_builtin_function
eyediagram(tx_Polar_Non_Return_Zero_signal,300,3);
title('Eye-Diagram for Polar Non Return to Zero');

% Transmit the signal
rx_Polar_Non_Return_Zero_signal= tx_Polar_Non_Return_Zero_signal;

% Sweep on the value of sigma and calculate BER for each value
sigma_values_2 = linspace(0, max(V_levels_2), 10);
ber_values_2 = zeros(1, length(sigma_values_2));
rx_bit_stream_2 = zeros (1,length(tx_bit_stream_2));
for i = 1:length(sigma_values_2)
    % Add noise to the received signal
    noise_2 = sigma_values_2(i) * randn(1, length(rx_Polar_Non_Return_Zero_signal));
    rx_Polar_Non_Return_Zero_signal_with_noise = rx_Polar_Non_Return_Zero_signal+ noise_2;

    for index = 50:100:length(rx_Polar_Non_Return_Zero_signal_with_noise)

        rx_bit_stream_2((index+50)/100) = (rx_Polar_Non_Return_Zero_signal_with_noise(index )> 0);

    end

    % Compare received bit stream with transmitted bit stream and count errors
    num_errors_2 = sum(rx_bit_stream_2 ~= tx_bit_stream_2);
    % Calculate BER for current value of sigma
    ber_values_2(i) = num_errors_2 / length(tx_bit_stream_2);
end

% Plot BER vs. sigma
figure (8)
hold all
semilogy(sigma_values_2, ber_values_2);
xlabel('Noise RMS Value (sigma)');
ylabel('Bit Error Rate (BER)');
title('BER vs. Noise RMS Value');
grid on;


%% Part_1_UPRZ


% Generating_random_stream_of_bits

tx_bit_stream_3 = randi([0, 1], 1, 10000); % 10,000 random bits

V_levels_3 = [0 1.2]; % Voltage levels for '0'and '1'
line_coded_stream_3 = V_levels_3(tx_bit_stream_3 + 1); % Assign voltage levels based on bit values

tx_Uni_polar_Return_Zero_signal = zeros(1, 1000000);
for Counter_3 = 1:length(line_coded_stream_3)
    Counter_1 = (Counter_3 - 1) * 100 + 1;
    Counter_2 = Counter_1 + 49;
    tx_Uni_polar_Return_Zero_signal(Counter_1:Counter_2) = line_coded_stream_3(Counter_3); % Set positive half of pulse
    Counter_1 = Counter_2 + 1;
    Counter_2 = Counter_1 + 49;
    tx_Uni_polar_Return_Zero_signal(Counter_1:Counter_2) = 0; % Set negative half of pulse
end


Counter_11=1;
Counter_22=100;
tx_bit_stream_33=zeros(1,1000000);
for Counter_33= 1:length(tx_bit_stream_3)
    for Counter_44 =Counter_11:Counter_22
        tx_bit_stream_33(Counter_44)=tx_bit_stream_3(Counter_33);
    end
    Counter_11=Counter_11+100;
    Counter_22=Counter_22+100;
end

% Time_Vector
Bits_duration =0.01;
t_3 = 0:Bits_duration:(length(tx_Uni_polar_Return_Zero_signal) * Bits_duration) - Bits_duration;


% Time-Domain Plot
figure (13)
subplot(2,1,1)
plot(t_3,tx_Uni_polar_Return_Zero_signal);
xlabel('Time');
ylabel('Voltage');
title('Uni-polar Return to Zero (RZ) Line Coding Example');
ylim([-0.5 1.5]);
xlim([0 5]);
grid on;

subplot(2,1,2)
plot(t_3,tx_bit_stream_33);
xlabel('Time');
ylabel('Value');
title('The Code');
xlim([0 5]);
ylim([-0.5 1.5]);
grid on;

%Frequency Domain Plot
Fs=1/Bits_duration;
T=10000;
df=1/T;
N=ceil(T/Bits_duration);
if (rem(N,2)==0)
    freq=(-0.5*Fs):df:(0.5*Fs)-df; %even function
else
    freq=-(0.5*Fs-0.5*df):df:(0.5*Fs-0.5*df); %odd function
end
freq_spectrum_3 = abs(fftshift(fft(tx_Uni_polar_Return_Zero_signal))).^2/N;
figure (9)
subplot(2,3,3)
plot(freq,freq_spectrum_3);
xlabel('Frequency');
ylabel('Magnitude squared');
title('Uni-polar Return to Zero');
grid on;



% Eye_diagram_plot_using_builtin_function
eyediagram(tx_Uni_polar_Return_Zero_signal,300,3);
title('Eye-Diagram for Uni-polar Return to Zero')

% Transmit the signal
rx_Uni_polar_Return_Zero_signal= tx_Uni_polar_Return_Zero_signal;

% Sweep on the value of sigma and calculate BER for each value
sigma_values_3 = linspace(0, max(V_levels_3), 10);
ber_values_3 = zeros(1, length(sigma_values_3));
rx_bit_stream_3 = zeros (1,length(tx_bit_stream_3));
for i = 1:length(sigma_values_3)
    % Add noise to the received signal
    noise_3 = sigma_values_3(i) * randn(1, length(rx_Uni_polar_Return_Zero_signal));
    rx_Uni_polar_Return_Zero_signal_with_noise = rx_Uni_polar_Return_Zero_signal+ noise_3;

    for index = 50:100:length(rx_Uni_polar_Return_Zero_signal_with_noise)

        rx_bit_stream_3((index+50)/100) = (rx_Uni_polar_Return_Zero_signal_with_noise(index )> 0.6);

    end

    % Compare received bit stream with transmitted bit stream and count errors
    num_errors_3 = sum(rx_bit_stream_3 ~= tx_bit_stream_3);
    % Calculate BER for current value of sigma
    ber_values_3(i) = num_errors_3 / length(tx_bit_stream_3);
end

% Plot BER vs. sigma
figure (8)
hold all
semilogy(sigma_values_3, ber_values_3);
xlabel('Noise RMS Value (sigma)');
ylabel('Bit Error Rate (BER)');
title('BER vs. Noise RMS Value');
grid on;


%% Part_1_BPRZ

% Generating_random_stream_of_bits
tx_bit_stream_4 = randi([0, 1], 1, 10000); % 10,000 random bits
V_levels_4 = [0 1.2]; % Voltage levels for '0'and '1'
line_coded_stream_4 = V_levels_4(tx_bit_stream_4 + 1); % Assign voltage levels based on bit values

tx_Bi_Polar_Return_Zero_Signal=zeros(1, 1000000);
Flagg=1.2;
for Counter_3 = 1:length(line_coded_stream_4)
    Counter_1 = (Counter_3 - 1) * 100 + 1;
    Counter_2 = Counter_1 + 49;
    if(line_coded_stream_4(Counter_3)>0)
        tx_Bi_Polar_Return_Zero_Signal(Counter_1:Counter_2) = Flagg;% Set positive half of pulse
        Flagg=-1*Flagg;
        Counter_1 = Counter_2 + 1;
        Counter_2 = Counter_1 + 49;
        tx_Bi_Polar_Return_Zero_Signal(Counter_1:Counter_2) = 0; % Set negative half of pulse
    end
end


Counter_11=1;
Counter_22=100;
tx_bit_stream_44=zeros(1,1000000);
for Counter_33= 1:length(tx_bit_stream_4)
    for Counter_44 =Counter_11:Counter_22
        tx_bit_stream_44(Counter_44)=tx_bit_stream_4(Counter_33);
    end
    Counter_11=Counter_11+100;
    Counter_22=Counter_22+100;
end

% Time_Vector
Bits_duration =0.01;
t_4 = 0:Bits_duration:(length(tx_Bi_Polar_Return_Zero_Signal) * Bits_duration) - Bits_duration;

% Time-Domain Plot
figure (14)
subplot(2,1,1)
plot(t_4,tx_Bi_Polar_Return_Zero_Signal);
xlabel('Time');
ylabel('Voltage');
title('Bi-polar return to zero (BPRZ) Line Coding Example');
ylim([-1.5 1.5]);
xlim([0 5]);
grid on;

subplot(2,1,2)
plot(t_4,tx_bit_stream_44);
xlabel('Time');
ylabel('Value');
title('The Code');
xlim([0 5]);
ylim([-0.5 1.5]);
grid on;

%Frequency Domain Plot
Fs=1/Bits_duration;
T=10000;
df=1/T;
N=ceil(T/Bits_duration);
if (rem(N,2)==0)
    freq=(-0.5*Fs):df:(0.5*Fs)-df; %even function
else
    freq=-(0.5*Fs-0.5*df):df:(0.5*Fs-0.5*df); %odd function
end
freq_spectrum_4 = abs(fftshift(fft(tx_Bi_Polar_Return_Zero_Signal))).^2/N;
figure (9)
subplot(2,3,4)
plot(freq,freq_spectrum_4);
xlabel('Frequency');
ylabel('Magnitude squared');
title('Bi-polar return to zero');
grid on;



% Eye_diagram_plot_using_builtin_function
eyediagram(tx_Bi_Polar_Return_Zero_Signal,300,3);
title('Eye-Diagram for Bi-polar return to zero');

% Transmit the signal
rx_Bi_Polar_Return_Zero_signal= tx_Bi_Polar_Return_Zero_Signal;

% Sweep on the value of sigma and calculate BER for each value
sigma_values_4 = linspace(0, max(V_levels_4), 10);
ber_values_4 = zeros(1, length(sigma_values_4));
errors_detected_Bi_Polar = zeros(1, length(sigma_values_4));
for i = 1:length(sigma_values_4)
    % Add noise to the received signal
    noise_4 = sigma_values_4(i) * randn(1, length(rx_Bi_Polar_Return_Zero_signal));
    rx_Bi_Polar_Return_Zero_signal_with_noise = rx_Bi_Polar_Return_Zero_signal+ noise_4;
    rx_bit_stream_4 = zeros (1,length(tx_bit_stream_4));
    rx_signed_bit_stream_4 = zeros (1,length(tx_bit_stream_4));

    for index = 50:100:length(rx_Bi_Polar_Return_Zero_signal_with_noise)

        rx_bit_stream_4((index+50)/100) = (abs (rx_Bi_Polar_Return_Zero_signal_with_noise(index ))> 0.6 );

        if ((rx_Bi_Polar_Return_Zero_signal_with_noise(index ))> 0.6 )
            rx_signed_bit_stream_4((index+50)/100) = 1;
        elseif ((rx_Bi_Polar_Return_Zero_signal_with_noise(index )) < -0.6 )
            rx_signed_bit_stream_4((index+50)/100) = - 1;
        end

    end
    pre=-1;
    % Error Detection in Bi-polar Return to Zero
    for b=1:length(rx_signed_bit_stream_4)
        if  rx_signed_bit_stream_4(b) ~= 0
            now =  rx_signed_bit_stream_4(b);
            if now == pre
                errors_detected_Bi_Polar(i) = errors_detected_Bi_Polar (i) + 1;
            end
            pre = now;
        end
    end

    % Compare received bit stream with transmitted bit stream and count errors
    num_errors_4 = sum(rx_bit_stream_4 ~= tx_bit_stream_4);
    % Calculate BER for current value of sigma
    ber_values_4(i) = num_errors_4 / length(tx_bit_stream_4);
    % errors_detected(i) = errors_detected(i) / length(tx_bit_stream_4);
end

% Plot BER vs. sigma
figure (8)
hold all
semilogy(sigma_values_4, ber_values_4);
xlabel('Noise RMS Value (sigma)');
ylabel('Bit Error Rate (BER)');
title('BER vs. Noise RMS Value');
grid on;

% Plot Errors vs. sigma
figure (10)
plot(sigma_values_4,ber_values_4*10000,sigma_values_4,errors_detected_Bi_Polar);
legend ('Real','Detected');
xlabel('Noise RMS Value (sigma)');
ylabel('Numbers of Errors');
title('Bi-polar return to zero (Real Errors vs. Errors Detected)');
grid on;


%% Part_1_Tx_MAN

% Generating_random_stream_of_bits

tx_bit_stream_5 = randi([0, 1], 1, 10000); % 10,000 random bits
V_levels_5 = [-1.2 1.2]; % Voltage levels for '0'and '1'
line_coded_stream_5 = V_levels_5(tx_bit_stream_5 + 1); % Assign voltage levels based on bit values

Counter_1=1;
Counter_2=100;
tx_MAN_Signal=zeros(1, 1000000);
for Counter_3 = 1:length(line_coded_stream_5)
    Counter_1 = (Counter_3 - 1) * 100 + 1;
    Counter_2 = Counter_1 + 49;
    tx_MAN_Signal(Counter_1:Counter_2) = line_coded_stream_5(Counter_3); % Set positive half of pulse
    Counter_1 = Counter_2 + 1;
    Counter_2 = Counter_1 + 49;
    tx_MAN_Signal(Counter_1:Counter_2) = -line_coded_stream_5(Counter_3); % Set negative half of pulse
end

Counter_11=1;
Counter_22=100;
tx_bit_stream_55=zeros(1,1000000);
for Counter_33= 1:length(tx_bit_stream_5)
    for Counter_44 =Counter_11:Counter_22
        tx_bit_stream_55(Counter_44)=tx_bit_stream_5(Counter_33);
    end
    Counter_11=Counter_11+100;
    Counter_22=Counter_22+100;
end

% Time_Vector
Bits_duration =0.01;
t_5 = 0:Bits_duration:(length(tx_MAN_Signal) * Bits_duration) - Bits_duration;

% Time-Domain Plot
figure (15)
subplot(2,1,1)
plot(t_5,tx_MAN_Signal);
xlabel('Time');
ylabel('Voltage');
title('Manchester Line Coding Example');
xlim([0 5]);
ylim([-1.5 1.5]);
grid on;


subplot(2,1,2)
plot(t_5,tx_bit_stream_55);
xlabel('Time');
ylabel('Value');
title('The Code');
xlim([0 5]);
ylim([-0.5 1.5]);
grid on;

%Frequency Domain Plot
Fs=1/Bits_duration;
T=10000;
df=1/T;
N=ceil(T/Bits_duration);
if (rem(N,2)==0)
    freq=(-0.5*Fs):df:(0.5*Fs)-df; %even function
else
    freq=-(0.5*Fs-0.5*df):df:(0.5*Fs-0.5*df); %odd function
end
freq_spectrum_5 = abs(fftshift(fft(tx_MAN_Signal))).^2/N;
figure (9)
subplot(2,3,5)
plot(freq,freq_spectrum_5);
xlabel('Frequency');
ylabel('Magnitude squared');
title('Manchester');
grid on;



% Eye_diagram_plot_using_builtin_function
eyediagram(tx_MAN_Signal,300,3);
title('Eye-Diagram for Manchester');

% Transmit the signal
rx_MAN_Signal= tx_MAN_Signal;

% Sweep on the value of sigma and calculate BER for each value
sigma_values_5 = linspace(0, max(V_levels_5), 10);
ber_values_5 = zeros(1, length(sigma_values_5));
rx_bit_stream_5 = zeros (1,length(tx_bit_stream_5));
for i = 1:length(sigma_values_5)
    % Add noise to the received signal
    noise_5 = sigma_values_5(i) * randn(1, length(rx_MAN_Signal));
    rx_MAN_Signal_with_noise = rx_MAN_Signal+ noise_5;

    for index = 50:100:length(rx_MAN_Signal_with_noise)

        rx_bit_stream_5((index+50)/100) = (rx_MAN_Signal_with_noise(index )> 0);

    end

    % Compare received bit stream with transmitted bit stream and count errors
    num_errors_5 = sum(rx_bit_stream_5 ~= tx_bit_stream_5);
    % Calculate BER for current value of sigma
    ber_values_5(i) = num_errors_5 / length(tx_bit_stream_5);
end

% Plot BER vs. sigma
figure (8)
hold all
semilogy(sigma_values_5, ber_values_5);
legend ("UPNRZ","PNRZ","UPRZ","BPRZ","MAN")
xlabel('Noise RMS Value (sigma)');
ylabel('Bit Error Rate (BER)');
title('BER vs. Noise RMS Value');
grid on;
