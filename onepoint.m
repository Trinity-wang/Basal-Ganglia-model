function Idbs=onepoint(Idbs,tempdbs,nlat,radius,j) 

position1=1;
position2=50;

wdbs = calculatedbsguass(nlat, tempdbs, radius);
    Idbs(position1:position2,position1:position2)=wdbs; 
%     figure(1)
%     mesh(Idbs);
%     pause(0.1)
end