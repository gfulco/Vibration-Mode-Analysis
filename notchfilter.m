function [ filter ] = notchfilter( position,order,type,width)

%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

fr = [-5120:0.1000:5120];


unit = fr>=0 | fr<0;

rect = fr>=(position-width) & fr<=(position+width);
rect_cc = fr>=-(position+width) & fr<=-(position-width);

ideal = unit - rect - rect_cc;

figure('name','ideal frequency domain filter ')
plot(fr, ideal);
title('Ideal Notch filter spectrum')
xlabel('f (Hz)')
ylabel('ideal(f)')

ideal_n = ifftshift(ifft(ideal,'symmetric'));
figure('name', 'ideal discrete domain filter')
plot([-5120:0.1:5120],ideal_n)
xlabel('n')
ylabel('ideal_n(n)')

fvtool(ideal_n)
window(1:length(ideal_n)) = 0;

triangle = triang(order);
rect = rectwin(order);
gauss = sqrt(8*pi)*normpdf(-order/20:0.1:order/20,0,2);
%gauss = sqrt(4*pi)*gauss/norm(gauss);

% figure('name', 'gauss')
% plot([1:order+1],gauss)
% xlabel('n')
% ylabel('window(n)')

if(strcmpi(type,'triangle'))
    window((end)/2:(end)/2+order) = triangle;
elseif(strcmpi(type,'gaussian'))
    window((end)/2:(end)/2+order) = gauss;
elseif(strcmpi(type,'rect'))
   window((end)/2:end/2+order) = rect;
end

figure('name', 'chosen window')
plot([-5120:0.1:5120],window)
ylabel('window(n)')


shifted = circshift(ideal_n,order/2);
N = numel(shifted);
ix = (1:N) - order/2;
tf = ix < 1 | ix > N;
shifted(tf) = 0 ;
shifted(1:51200+1) = 0;
shifted(51200+order+1:end) = 0;

figure('name', 'shift')
hold on
plot([-5120:0.1:5120],shifted)
plot([-5120:0.1:5120],window)


windowed = shifted.*window;
figure('name', 'windowed filter')
plot([-5120:0.1:5120],windowed)


filter = fftshift(fft(windowed,length(fr)));
filter = 51200*filter;
 
scaled = abs(filter/length(fr));
temp = scaled(1:length(fr)/2+1); 
temp(2:end-1) = 2*temp(2:end-1);
figure('name','frequency real filter')
plot([0:0.1:5120],temp) 
title('Single-Sided Amplitude Spectrum of Filter')
xlabel('f (Hz)')
ylabel('|Filter(f)|')





%z = ztransf(windowed);

end

