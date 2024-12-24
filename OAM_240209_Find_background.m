
% funtion to select the background of leaves. 

function [bck, bck_sampler]=OAM_240209_Find_background(I)


%numbM=round(size(f_names,2));
cord=zeros(20,2); % I=c8;
figure(1);imagesc(I);
% f1.WindowState = 'maximized';
title([ 'select background areas ']);%'pos =' num2str(pos)

button1=0;
A=0; 
while button1~=1000
[X,Y,button]=ginput(1);

X=round(X); Y=round(Y);

if button(1) == 1
cord(A+1,1)=X(1);
cord(A+1,2)=Y(1);
A=A+1;
elseif button(1) ==2 
button1=1000;
end 
end

close

numbM=sum(logical(cord(:,1)~=0));
med=zeros();
bck_sampler=zeros();
interv=1:2:numbM;

for ig =interv%size(f_names,2) % [1:343 345:790]% itt=1
IX=double(I(cord(ig,2):cord(ig+1,2),cord(ig,1):cord(ig+1,1)));% 
%figure;imagesc(IX)
%pause
med(1,interv==ig)=median(IX(:),"omitmissing");% <=

if ig==1
bck_sampler=reshape(IX,[numel(IX),1]);
else
bck_sampler=[bck_sampler;reshape(IX,[numel(IX),1])];
end

close


end

bck=mean(med);

