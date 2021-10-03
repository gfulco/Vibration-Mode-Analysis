function [ y ] = movav( sig, f , fc, l )
% Apply a moving average filter to sig
% f frequency axis vector
% l  window for the moving avg

y = sig*0;

for i = 1:length(sig)
    temp = 0;
    if i<=l
        for k = 1:i   
            temp = temp + sig(k);
         end
         y(i) = temp/i;
    elseif (i+l)>length(sig)
       for k = (i-l):length(sig)
           temp = temp+sig(k);
       end
        y(i) = temp/(length(sig)-i+l+1);
    elseif i > l
        for k = (i-l):(i+l)
            temp = temp+sig(k);
        end
        y(i) = temp/(2*l+1);
    end
end

 figure('Name','Moving Avg filter applied to signal')
 plot(f(1:fc*10),y(1:fc*10));

end

