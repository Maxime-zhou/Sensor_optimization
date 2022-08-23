% checking the orientation vector L to ensure the satisfy the constrains：
% 1. fiber path can't go back (overlapped)
% 2. the turinng angle of each segment should not less than 90°

function c=Check_L(L)
c = 1;
for i = 1:length(L)-1
%     if abs(L(i+1)-L(i))==3 || abs(L(i+1)-L(i))==4 || abs(L(i+1)-L(i))==5
    if abs(L(i+1)-L(i))==3 || abs(L(i+1)-L(i)) ==5  
    c = 0; 
    end
end