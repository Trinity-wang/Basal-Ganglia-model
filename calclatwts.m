function w = calclatwts(nlat, sig, rad)
%nlat is the size of the window. It must be an odd number.
%sig: max strength of the connections
%rad: radius of neighborhood

w = zeros(nlat, nlat);
ic = (nlat+1)/2;
jc = (nlat+1)/2;


for i = 1:nlat,
   for j = 1:nlat,
        dis = (i-ic)*(i-ic) + (j-jc)*(j-jc);
        w(i, j) = sig*exp(-dis/(rad*rad)) ; 
        if(i==j)
            w(i,j)=0;
        end
   end
end



