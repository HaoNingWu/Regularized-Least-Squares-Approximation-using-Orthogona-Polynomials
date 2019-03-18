function h = Filter(L,l)
if l/L <=0.5
    h = 1;
elseif 1/L <=1
    h = sin(pi*l/L)^2;
else
    h = 0;
end
% if l/L <= 0.25
%     h = 1.5;
% elseif l/L <= 0.75
%     h = .2;
% elseif l/L<=1
%     h = 1.5;
% else
%     h = 0;
% end
end
    