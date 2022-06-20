%% ROFdenoise
%
% Copyright (c) 2022, Carl Löndahl & Philippe Magiera
%
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are met:
%
% Redistributions of source code must retain the above copyright notice, 
% this list of conditions and the following disclaimer.
% Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.

% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
% I MPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE 
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF 
% THE POSSIBILITY OF SUCH DAMAGE
%
%  Description:
%
%  This denoising method is based on total-variation, originally proposed by
%  Rudin, Osher and Fatemi. In this particular case fixed point iteration
%  is utilized.
%
%  For the included image, a fairly good result is obtained by using a
%  theta value around 12-16. A possible addition would be to analyze the
%  residual with an entropy function and add back areas that have a lower
%  entropy, i.e. there are some correlation between the surrounding pixels.
% 


function A = ROFdenoise(Image, Theta)
[Image_h Image_w] = size(Image); 
g = 1; dt = 1/4; nbrOfIterations = 5;
Image = double(Image);
p = zeros(Image_h,Image_w,2);
d = zeros(Image_h,Image_w,2);
div_p = zeros(Image_h,Image_w);
for i = 1:nbrOfIterations
    for x = 1:Image_w
        for y = 2:Image_h-1
            div_p(y,x) = p(y,x,1) - p(y-1,x,1);
        end
    end
    for x = 2:Image_w-1
        for y = 1:Image_h
            div_p(y,x) = div_p(y,x) + p(y,x,2) - p(y,x-1,2);
        end
    end
    
    % Handle boundaries
    div_p(:,1) = p(:,1,2);
    div_p(:,Image_w) = -p(:,Image_w-1,2);
    div_p(1,:) = p(1,:,1);
    div_p(Image_h,:) = -p(Image_h-1,:,1);
    % Update u
    u = Image-Theta*div_p;
    % Calculate forward derivatives
    du(:,:,2) = u(:,[2:Image_w, Image_w])-u;
    du(:,:,1) = u([2:Image_h, Image_h],:)-u;
    % Iterate
    d(:,:,1) = (1+(dt/Theta/g).*abs(sqrt(du(:,:,1).^2+du(:,:,2).^2)));
    d(:,:,2) = (1+(dt/Theta/g).*abs(sqrt(du(:,:,1).^2+du(:,:,2).^2)));
    p = (p-(dt/Theta).*du)./d;
    
end
A = u;