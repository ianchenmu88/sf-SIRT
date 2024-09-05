
clear all
close all
load recPixels

s = radon(recPixels,[25:154]);


theta=[25:154];

C = fft(s); C1 = C;
cf = round(size(C, 1)/2);
SS = sum(abs(C).^2, 2);

%SS1 = SS(1:cf) + SS(cf:end); 
cc  = 5;
SS1 = SS(1:cf-cc);
SS1 = SS1 + SS(end:-1:cf+cc);
%SS1 = SS(1:cf);
%SS1 = SS1 + SS(end:-1:cf);

[SSs,id] = sort(SS1, 'descend'); 

d       = length(SSs);
n       = size(s, 1);

% gMDL
RSS = (sum(SSs) - cumsum(SSs'))/length(theta);
S1  = RSS ./ (n-2*(1:d));
F1  = (cumsum(SSs')/length(theta)) ./ (2*(1:d).*S1);
gMDL = log(S1) + 0.5*((1:d)/n).*log(F1);
[~, d] = min(gMDL(1:round(d/2)));


C1(cc-1+id(d+1:end), :) = 0; 
C1(size(C, 1)-cc-id(d+1:end), :) = 0;
Rhat = ifft(C1);
Rhat = real(Rhat);


At = iradon(Rhat,theta,'linear', 'cosine', 1,345); %reconstruct noisy alien

n = 10;%iterations
Fk = At;%Matrix Fk is our solution at the k-th step, now it is our initial guess
for  k=1:n
%    imshow(Fk, []);    
    R=   radon(Fk,theta);
    C = fft(R); 
    cf = round(size(C, 1)/2);
    SS = sum(abs(C).^2, 2);

    SS1 = SS(1:cf) + SS(cf:end); 
    SS1 = SS(1:cf-cc);
    SS1 = SS1 + SS(end:-1:cf+cc);
    %SS1 = SS(1:cf);
    %SS1 = SS1 + SS(end:-1:cf);

    [SSs,id] = sort(SS1, 'descend'); 

    d       = length(SSs);
    n       = size(R, 1);
    % gMDL
    RSS = (sum(SSs) - cumsum(SSs'))/length(theta);
    S1  = RSS ./ (n-2*(1:d));
    F1  = (cumsum(SSs')/length(theta)) ./ (2*(1:d).*S1);
    gMDL = log(S1) + 0.5*((1:d)/n).*log(F1);
    [~, d] = min(gMDL(1:round(d/2)));


    C(cc-1+id(d+1:end), :) = 0; 
    C(size(C, 1)-cc-id(d+1:end), :) = 0;
    Rhat = ifft(C);
    Rhat = real(Rhat);

    t = iradon(Rhat,theta, 'linear', 'cosine', 1,345);% reconstruct alien using Fk unfiltered sinogram

    Fk = Fk + At - t;
 
end

imshow(Fk,[])
