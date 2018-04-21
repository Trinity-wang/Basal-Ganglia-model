function [action,i]=selection(I)
z(:,1)=zeros(2,1);
threshold=0.15;
dt=0.01;
[a,b]=size(I);
b;
for i=2:b
    z(:,i)=z(:,i-1).*(1-dt)+dt.*I(:,i-1);
    
    Z=z(:,i);
    
    II= I(:,i-1);
    
    %pause();
    [M, index]=(max(z(:,i)));
    if M>threshold
%         i
      
        action=index;
    break;
    else 
        action =0;
    end
end
action;
end    