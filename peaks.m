function [ min_loc, max_loc ] = peaks( y, f, fc, L, k, p )
% Find peaks of a transfer function
% y amplitude vector
% f frequency axis vector
% L  minimum peak distance
% fc  Max frequency we want to analyze
% k is the width of the peak i want to keep
% p is the minimum increase in the next sample to be considered as a peak
j = 1;
%p = 0.0024
fc_index = find(f<fc,1,'last');

slope = sign(diff(y));
sign_switch(1:fc_index) = 0;

for i = 1:fc_index
        if(slope(i) > slope(i+1))
        sign_switch(i)=y(i+1);
        else
        sign_switch(i)=0;
        end
end

while j > 0
    temp = max(sign_switch);
    m = find(y == temp,1,'first');
    if(temp == 0)
        j=0;
    else
            maximum(1,j) = m-0.1*k;
            maximum(2,j) = m+0.2*k;
   
        
        if(maximum(1,j) > L && maximum(1,j)<(fc_index-L))
            for i = (maximum(1,j)-L):(maximum(1,j)+L)
                sign_switch(i) = 0;
            end
        elseif (maximum(1,j)<L)
            for i = 1:(maximum(1,j)+L)
                sign_switch(i) = 0;
            end
        elseif (maximum(1,j)>fc_index-L)
            for i = (maximum(1,j)-L):fc_index
                sign_switch(i) = 0;
            end
        end
        j = j+1;
    end
end

for b = 1:length(maximum)
    
        if(maximum(1,b)-k > 0 && maximum(2,b)+k <= length(y))
             m = mean(y(maximum(1,b)-k:maximum(2,b)+k));
             g = y(maximum(1,b)+0.1*k);
             if(m > g)
                maximum(1,b) = 0;
                maximum(2,b) = 0;
             end
        elseif (maximum(2,b)+k > length(y))
            m = mean(y(maximum(1,b)-k:length(y)));
            g = y(maximum(1,b)+0.1*k);
            if(m > g)
                maximum(1,b) = 0;
                maximum(2,b) = 0;
             end
        elseif (maximum(1,b)-k < 1)
            m = mean(y(1:maximum(2,b)+k));
            g = y(maximum(1,b)+0.1*k);
            if(m > g)
                maximum(1,b) = 0;
                maximum(2,b) = 0;
             end
        end
       
end
maximum( :, ~any(maximum,1) ) = []; %delete 0 from matrix by column
temp(1,1) = 1;
while temp(1,1)>0
    temp(1,1) = 0;
    for ii = 1:length(maximum)-1
        if(maximum(1,ii) > maximum(1,ii+1))
            temp(1,1) = maximum(1,ii+1);
            temp(2,1) = maximum(2,ii+1);
            maximum(1,ii+1) = maximum(1,ii);
            maximum(2,ii+1) = maximum(2,ii);
            maximum(1,ii) = temp(1,1);
            maximum(2,ii) = temp(2,1);    
        end
    end
end

min_loc = maximum(1,:);
max_loc = maximum(2,:);
end