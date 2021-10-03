clear all
close all
clc

%% Data Loading


data_1=dlmread('hammer_data_1.txt');
data_2=dlmread('hammer_data_2.txt');
data_3=dlmread('hammer_data_3.txt');
data_4=dlmread('hammer_data_4.txt');

time = data_1(:,1);
%forces
force(1,:) = data_1(:,2);
force(2,:) = data_2(:,2);
force(3,:) = data_3(:,2);
force(4,:) = data_4(:,2);
%accelerometer 1

resp1(1,:) = data_1(:,3);
resp1(2,:) = data_2(:,3);
resp1(3,:) = data_3(:,3);
resp1(4,:) = data_4(:,3);
%accelerometer 2
resp2(1,:) = data_1(:,4);
resp2(2,:) = data_2(:,4);
resp2(3,:) = data_3(:,4);
resp2(4,:) = data_4(:,4);
%accelerometer 3
resp3(1,:) = data_1(:,5);
resp3(2,:) = data_2(:,5);
resp3(3,:) = data_3(:,5);
resp3(4,:) = data_4(:,5);

%% Fourier Transforms, Frequency and Impulsive Responses

L = length(time);                 %Time axis length
Fs = 10240;                       %Sampling Rate
f = (Fs*(0:(L/2))/L);             %Creating axis of frequency in Hz

%Forces

for j = 1:4
fourier_force(j,:) = fft(force(j,:));
end

y = dirac(time);
idx = y == Inf;
y(idx) = 0;
idx = 32646;
y(idx) = 45.33;

figure('Name','Real impulses')
subplot(4,1,1)
plot(time,force(1,:));
axis([0 3.5 -0.1 50])
subplot(4,1,2)
%plot(time,y);
plot(time,force(2,:));
axis([0 3.5 -0.1 50])
subplot(4,1,3)
plot(time,force(3,:));
axis([0 3.5 -0.1 50])
subplot(4,1,4)
plot(time,force(4,:));
axis([0 3.5 -0.1 50])


% %Accelerometer1-
fourier_acc1(1:4,1:length(resp1)) = 0;
freq_resp1(1:4,1:length(fourier_acc1)) = 0;
imp_resp1(1:4,1:length(freq_resp1)) = 0;
mean_fr1(1:length(freq_resp1)) = 0;

for j = 1:4
 fourier_acc1(j,:) = fft(resp1(j,:));
 freq_resp1(j,:) = fourier_acc1(j,:)./fourier_force(j,:);
 imp_resp1(j,:) = ifft(freq_resp1(j,:));
 freq_resp1(j,:) = abs(freq_resp1(j,:));
 mean_fr1 = mean_fr1 + freq_resp1(j,:)./4;
 
end


fourier_acc2(1:4,1:length(resp2)) = 0;
freq_resp2(1:4,1:length(fourier_acc2)) = 0;
imp_resp2(1:4,1:length(freq_resp2)) = 0;
mean_fr2(1:length(freq_resp2)) = 0;
for j = 1:4
 fourier_acc2(j,:) = fft(resp2(j,:));
 freq_resp2(j,:) = fourier_acc2(j,:)./fourier_force(j,:);
 imp_resp2(j,:) = ifft(freq_resp2(j,:));
 freq_resp2(j,:) = abs(freq_resp2(j,:));
 mean_fr2 = mean_fr2 + freq_resp2(j,:)./4;
end


fourier_acc3(1:4,1:length(resp3)) = 0;
freq_resp3(1:4,1:length(fourier_acc3)) = 0;
imp_resp3(1:4,1:length(freq_resp3)) = 0;
mean_fr3(1:length(freq_resp3)) = 0;

for j = 1:4
 fourier_acc3(j,:) = fft(resp3(j,:));
 freq_resp3(j,:) = fourier_acc3(j,:)./fourier_force(j,:);
 imp_resp3(j,:) = ifft(freq_resp3(j,:));
 freq_resp3(j,:) = abs(freq_resp3(j,:));
 mean_fr3 = mean_fr3 + freq_resp3(j,:)./4;
end



%% Plotting frequency and impulsive responses

figure('Name','Impulse responses for accelerometer 1')
for j= 1:4
subplot(4,1,j)
plot(time,imp_resp1(j,:));
axis([0 0.05 -4 4])
end
figure('Name','Impulse responses for accelerometer 2')
for j= 1:4
subplot(4,1,j)
plot(time,imp_resp2(j,:));
axis([0 0.05 -4 4])
end

figure('Name','Impulse responses for accelerometer 3')
for j= 1:4
subplot(4,1,j)
plot(time,imp_resp3(j,:));
axis([0 0.05 -4 4])
end

meanimp1 = imp_resp1(1,:)+imp_resp1(2,:)+imp_resp1(3,:)+imp_resp1(4,:);
meanimp1 = meanimp1/4;
meanimp2 = imp_resp2(1,:)+imp_resp2(2,:)+imp_resp2(3,:)+imp_resp2(4,:);
meanimp2 = meanimp2/4;
meanimp3 = imp_resp3(1,:)+imp_resp3(2,:)+imp_resp3(3,:)+imp_resp3(4,:);
meanimp3 = meanimp3/4;
figure('Name','Mean Impulse response for accelerometers')
subplot(3,1,1);
plot(time,meanimp1);
axis([0 .1 -5 5]);
subplot(3,1,2);
plot(time,meanimp2);
axis([0 .1 -5 5]);
subplot(3,1,3);
plot(time,meanimp3);
axis([0 .1 -5 5]);

% %Single sided spectrum for mean responses

fr_sss(1:3,1:L/2+1) = 0;

fr_sss(1,:) = mean_fr1(1:L/2+1);
fr_sss(1,2:end-1) = 2*fr_sss(1,2:end-1);

fr_sss(2,:) = mean_fr2(1:L/2+1);
fr_sss(2,2:end-1) = 2*fr_sss(2,2:end-1);

fr_sss(3,:) = mean_fr3(1:L/2+1);
fr_sss(3,2:end-1) = 2*fr_sss(3,2:end-1);

%plotting the 3 spectrum
figure ('Name','Accelerometer 1 mean frequency response');
subplot(3,1,1);

plot(f,fr_sss(1,:));
axis([0 2500 0 500]);
subplot(3,1,2);
%figure ('Name','Accelerometer 2 mean frequency response');
plot(f,fr_sss(2,:));
axis([0 2500 0 500]);
subplot(3,1,3);
%figure ('Name','Accelerometer 3 mean frequency response');
plot(f,fr_sss(3,:));
axis([0 2500 0 500]);

%% finding the peaks


fc = 2000; %given, finding peaks until 2k Hz
avg_window = 260;
peak_range = 1270;
peak_width = 3200;

filtered(1:12,1:length(freq_resp1)) = 0;

for j = 1:4
    filtered(j,:) = movav(freq_resp1(j,:),f,fc,avg_window);
    filtered(j+4,:) = movav(freq_resp2(j,:),f,fc,avg_window);
    filtered(j+8,:) = movav(freq_resp3(j,:),f,fc,avg_window);
    
end    


figure('Name','Frequency response / filtered comparison');
hold on
plot(1:0.1:2000,freq_resp1(1,1:19991))
plot(1:0.1:2000,filtered(1,1:19991))
hold off

% 
for j = 1:4
    [min_loc(j,:),max_loc(j,:)] = peaks(filtered(j,:),f,fc,peak_width,peak_range);
end
[min_loc(5,:),max_loc(5,:)] = peaks(filtered(5,:),f,fc,3000,1800);
[min_loc(6,:),max_loc(6,:)] = peaks(filtered(6,:),f,fc,peak_width,2000);
[min_loc(7,:),max_loc(7,:)] = peaks(filtered(7,:),f,fc,peak_width,2000);
[min_loc(8,:),max_loc(8,:)] = peaks(filtered(8,:),f,fc,peak_width,2000);

for j = 1:4
    [min_loc(j+8,:),max_loc(j+8,:)] = peaks(filtered(j+8,:),f,fc,peak_width,peak_range);
end 

[row,col] = size(min_loc);

amplitude(1:row,1:col) = 0;
position(1:row,1:col) = 0;

for i = 1:4
    for j = 1:4
        amplitude(j,i) = max(freq_resp1(j,min_loc(j,i):max_loc(j,i)));
        
        amplitude(j+4,i) = max(freq_resp2(j,min_loc(j+4,i):max_loc(j+4,i)));
        amplitude(j+8,i) = max(freq_resp3(j,min_loc(j+8,i):max_loc(j+8,i)));
         
        position(j,i) = min_loc(j,i) + find(freq_resp1(j,min_loc(j,i):max_loc(j,i)) == amplitude(j,i),1,'first');
        position(j+4,i) = min_loc(j+4,i) + find(freq_resp2(j,min_loc(j+4,i):max_loc(j+4,i)) == amplitude(j+4,i),1,'first');
        position(j+8,i) = min_loc(j+8,i) + find(freq_resp3(j,min_loc(j+8,i):max_loc(j+8,i)) == amplitude(j+8,i),1,'first');
    
        position(j,i) = f(position(j,i));
        position(j+4,i) = f(position(j+4,i));
        position(j+8,i) = f(position(j+8,i));
        
        figure('Name','aaa')
        hold on
        plot(f(min_loc(j,i):max_loc(j,i)),freq_resp1(j,min_loc(j,i):max_loc(j,i)))
        stem(position(j,i),amplitude(j,i));
        hold off
    end
end

freq_resp_ss(1,:) = freq_resp1(1,1:L/2+1);
freq_resp_ss(1,2:end-1) = freq_resp_ss(1,2:end-1);

freq_resp_ss(2,:) = freq_resp1(2,1:L/2+1);
freq_resp_ss(2,2:end-1) = freq_resp_ss(1,2:end-1);

freq_resp_ss(3,:) = freq_resp1(3,1:L/2+1);
freq_resp_ss(3,2:end-1) = freq_resp_ss(1,2:end-1);

freq_resp_ss(4,:) = freq_resp1(4,1:L/2+1);
freq_resp_ss(4,2:end-1) = freq_resp_ss(1,2:end-1);

% freq_resp_ss(2,:) = freq_resp2(1,1:L/2+1);
% freq_resp_ss(2,2:end-1) = freq_resp_ss(2,2:end-1);
% 
% freq_resp_ss(3,:) = freq_resp3(1,1:L/2+1);
% freq_resp_ss(3,2:end-1) = freq_resp_ss(3,2:end-1);


figure('Name','Peaks found')
subplot(4,1,1)
hold on
plot(f,freq_resp_ss(1,:));
stem(position(1,1:4),amplitude(1,1:4))
axis([0 2100 -0.1 500])
hold off

subplot(4,1,2)
hold on
plot(f,freq_resp_ss(2,:));
stem(position(2,1:4),amplitude(2,1:4))
axis([0 2100 -0.1 500])
hold off

subplot(4,1,3)
hold on
plot(f,freq_resp_ss(3,:));
stem(position(3,1:4),amplitude(3,1:4))
axis([0 2100 -0.1 500])
hold off

subplot(4,1,4)
hold on
plot(f,freq_resp_ss(4,:));
stem(position(4,1:4),amplitude(4,1:4))
axis([0 2100 -0.1 500])
hold off
% subplot(3,1,2)
% hold on
% plot(f,freq_resp_ss(2,:));
% stem(position(5,1:4),amplitude(5,1:4))
% axis([0 2100 -0.1 500])
% hold off
% subplot(3,1,3)
% hold on
% plot(f,freq_resp_ss(3,:));
% stem(position(9,1:4),amplitude(9,1:4))
% axis([0 2100 -0.1 500]) 
% hold off




%% mean value and std deviation of the peaks

mean_pos(1:col) = 0;
mean_amp(1:col) = 0;
std_pos(1:col) = 0;
std_amp(1:col) = 0;

for j = 1:col
    for i = 1:row
        mean_pos(j) = mean_pos(j)+position(i,j)/row;
        mean_amp(j) = mean_amp(j)+amplitude(i,j)/row;
    end
end

for j = 1:col
    for i = 1:row 
         std_pos(j) = std_pos(j) + (position(i,j)-mean_pos(j))^2;
         std_amp(j) = std_amp(j) + (amplitude(i,j)-mean_amp(j))^2;
    end
    std_pos(j) = sqrt(std_pos(j)/(row-1));
    std_amp(j) = sqrt(std_amp(j)/(row-1));
end

fre_mean =( fr_sss(1,:)+fr_sss(2,:)+fr_sss(3,:))./3;
figure('Name','Mean position of peaks')
hold on
plot(f,fre_mean(:));
stem(mean_pos(1),2*mean_amp(1))
stem(mean_pos(2),2*mean_amp(2))
stem(mean_pos(3),2*mean_amp(3))
stem(mean_pos(4),2*mean_amp(4))
% stem(mean_pos(3)-std_pos(3),2*(mean_amp(3)-std_amp(3)),'g')
% stem(mean_pos(3)-std_pos(3),2*(mean_amp(3)+std_amp(3)),'g')
% stem(mean_pos(3)+std_pos(3),2*(mean_amp(3)+std_amp(3)),'r')
% stem(mean_pos(3)+std_pos(3),2*(mean_amp(3)-std_amp(3)),'r')
axis([0 2100 -0.1 500]) 
hold off

close all
clc
%%  notch filter design

filter = notchfilter(mean_pos(1),500,'gaussian',mean_pos(1)/2);
% tri = freq_resp1(1,:).*filter(1:102400);
% 
% fil(1,:) = abs(tri(1:L/2+1)/L);
% fil(1,2:end-1) = 2*fil(1,2:end-1);

scaled = abs(filter/L);
temp = scaled(1:L/2+1); 
temp(2:end-1) = 2*temp(2:end-1);

tri = fr_sss(1,:).*temp(1,:);


% 
% res(1,:) = abs(freq_resp1(1,1:L/2+1)/L);
% res(1,2:end-1) = 2*res(1,2:end-1);
figure ('Name','Filtered');
hold on
plot(f,tri(1,:));
plot(f,fr_sss(1,:));
axis([0 2000 0 800]);
hold off

%%



