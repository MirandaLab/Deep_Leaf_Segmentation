


% load toydata:

I=imread('Colletotrichum_orbiculare_infected_leaves.jpg');

% seect are to be processed
c0=double(I(490:2650,1000:3230,:));
% offset the image to avoid negative pixels after background substraction
c0=c0+100;
c0=(c0-min(c0(:)))./(max(c0(:))-min(c0(:)));
% visualize:
%  imshow([c0]); 

c1=double(c0(:,:,1));
c1=(c1-min(c1(:)))./(max(c1(:))-min(c1(:))); % imshow([c1]); imagesc(c1,[0.1 0.2])

c2=double(c0(:,:,2));
c2=(c2-min(c2(:)))./(max(c2(:))-min(c2(:)));% imshow([c2]);imagesc(c2,[0.1 0.2])

c3=double(c0(:,:,3));
c3=(c3-min(c3(:)))./(max(c3(:))-min(c3(:)));% imshow([c3]);imagesc(c3,[0.1 0.2])

%visualize
%  imagesc([c1 c2 c3])

%% select background pixels from each spectral channels
% click the upper left and lower right corner delimiting a squared area in
% the background

A2=zeros(size(c0(:,:,1)));
for iy=1:3
[bck, bck_sampler]=OAM_240209_Find_background(c0(:,:,iy));
A1=c0(:,:,iy)-2*bck;
imshow(imfill(bwmorph(double(imbinarize(A1)),'thicken',5),"holes"))
A2=A2+A1;
%pause
end

%% remove pixels of lower value than the background
A3=A2;
A3(A3<=0.5*median(A3(:)))=0; % imagesc(A3)

% imshow(A3)

%% binarize the pixels for the preliminary mask
A4=A3;
A3(A3~=0)=1;
% imshow([A3])

%% eliminate pixels selected outside the leaves
A5=bwareafilt(logical(A3),[1 10000]);% imshow([A5])
A6=(double(~A5).*A3);% imshow([A6])
A7=bwlabel(imbinarize(A6));% 
%imagesc([A7]); % inshow doesn't display the masks

%% separate leaves into individual images
% first row: binary mask
% second row original image
% third row flatfield corrected image

numbf=max(A7(:));
leavs=cell(3,numbf);
for io=1:max(A7(:))
    A8=(A7==io); 
    A9=imfill(A8,"holes");
    [x2,y2]=get_wind_coord1(A9,1);
    leavs{1,io}=A9(y2,x2);% mask % imagesc(leavs{1,io})
    leavs{2,io}=c0(y2,x2,:);%  mask % imagesc(leavs{2,io})
    leavs{3,io}=imflatfield(leavs{2,io},10,leavs{1,io}); % imagesc(leavs{3,io}) 
end

for io1=1:max(A7(:))
imagesc([cat(3,leavs{1,io1},leavs{1,io1},leavs{1,io1}) leavs{2,io1} leavs{3,io1}] )
pause
end


%% quantify the number of preliminary identified as infected vs the total number
% of pixels identified as part of the leaf
% using the original the flatfield corrected, or a mixture of both.

  quanti0=zeros(2,2);
  interv=[3 4 6 7 10 11 15 14 19 17; 1 2 5 6 9 8 12 13 16 18]; % first control, second treatment
 quanti1=cell(2,1);
 sae=1;


for aa=1:2

  quanti=zeros(numel(interv(aa,:)),2); 
for io2=interv(aa,:)%  io2 =3
 
    if io2>19
        continue
    end

A1B=rgb2gray(imerode(leavs{1,io2},ones(9)).*leavs{2,io2});% figure;imagesc(A1B)
A2B=rgb2gray(imerode(leavs{1,io2},ones(9)).*leavs{3,io2});% figure;imagesc(A2B)

% A1=leavs{1,io2}.*leavs{2,io2};% imagesc(A1)
% A2=leavs{1,io2}.*leavs{3,io2};% imagesc(A2)
% A1a=A1(:,:,2);
% A2a=A2(:,:,2);

med1=median(A1B(A1B~=0));
med2=median(A2B(A2B~=0));

A3B=imopen(imbinarize(A1B-1*med1),ones(4)); % figure;imshow(A3B)
A4B=imopen(imbinarize(A2B-1*med2),ones(4)); % figure;imshow(A4B)

% A3a=A3.*leavs{2,io2};%figure;imshow(A3a)
A3b=cat(3,(0.9.*A3B+leavs{2,io2}(:,:,1)), leavs{2,io2}(:,:,2), leavs{2,io2}(:,:,3));%  figure(1);imshow(A3b)
inf_a=numel(find(A3B~=0))./numel(find(leavs{1,io2}~=0));

% A4a=A4.*leavs{3,io2};figure(2);%imshow(A4a); imshow(leavs{1,io2})
A4b=cat(3,(0.9.*A4B+leavs{3,io2}(:,:,1)), leavs{3,io2}(:,:,2), leavs{3,io2}(:,:,3));%  figure(2);imshow(A4b)
inf_b=numel(find(A4B~=0))./numel(find(leavs{1,io2}~=0));


A5b=imbinarize(A4B+A3B); 

leavs{4,io2}=A5b;% imagesc(leavs{4,io2})
  

A5c=cat(3,(0.9.*A5b+leavs{3,io2}(:,:,1)), leavs{3,io2}(:,:,2), leavs{3,io2}(:,:,3));%figure(2);imshow(A4b)
inf_c=numel(find(A5b~=0))./numel(find(leavs{1,io2}~=0));


io3=find(interv(aa,:)==io2);
quanti(io3,1)=inf_a*100;
quanti(io3,2)=inf_b*100;
quanti(io3,3)=inf_c*100;

f1=figure(1);imshow([leavs{2,io2} A3b A4b A5c]);
text(size(leavs{2,io2},2),10,['Affected area = ' num2str(round(quanti(io3,1))) '%'],"FontSize",11,"FontWeight","bold","Color",[1 1 1]);
text(size(leavs{2,io2},2)*2,10,['Affected area = ' num2str(round(quanti(io3,2))) '%'],"FontSize",11,"FontWeight","bold","Color",[1 1 1]);
text(size(leavs{2,io2},2)*3,10,['Affected area = ' num2str(round(quanti(io3,3))) '%'],"FontSize",11,"FontWeight","bold","Color",[1 1 1]);

if sae==1
saveas(f1,['Leave_no' num2str(io2)],'pdf')
end

% pause (0.05)

% imagesc([cat(3,leavs{1,io1},leavs{1,io1},leavs{1,io1}) leavs{2,io1} leavs{3,io1}] )
% pause
end

quanti(quanti==0)=nan;

quanti0(aa,1)=mean(quanti(1:max(find(quanti(:,1)~=0)),1),"omitmissing");
quanti0(aa,2)=std(quanti(1:max(find(quanti(:,1)~=0)),1),"omitmissing");
quanti0(aa,3)=mean(quanti(1:max(find(quanti(:,1)~=0)),2),"omitmissing");
quanti0(aa,4)=std(quanti(1:max(find(quanti(:,1)~=0)),2),"omitmissing");
quanti0(aa,5)=mean(quanti(1:max(find(quanti(:,1)~=0)),3),"omitmissing");
quanti0(aa,6)=std(quanti(1:max(find(quanti(:,1)~=0)),3),"omitmissing");

quanti1{aa,1}=quanti;

end



% calculate the significance of the percentages of affected area
pval=zeros(1,3);
for ee=1:3
[a,b]=kstest2(quanti1{1,1}(:,ee),quanti1{2,1}(:,ee));
pval(1,ee)=b;
end

% visualized both images and masks
% for io1=1:max(A7(:))
% imagesc([cat(3,leavs{1,io1},leavs{1,io1},leavs{1,io1}) leavs{2,io1} leavs{3,io1} cat(3,leavs{4,io1},leavs{4,io1},leavs{4,io1})] )
% pause
% end

% visualize
% imagesc([leavs{4,io1},leavs{4,io1},leavs{4,io1}])


bar([quanti0(:,1)' quanti0(:,3)' quanti0(:,5)'])
hold on
errorbar(1:6,[quanti0(:,1)' quanti0(:,3)' quanti0(:,5)'],[quanti0(:,2)' quanti0(:,4)' quanti0(:,6)'],"o")
ylim([0 25])
legend(num2str([pval(1,1) pval(1,2) pval(1,3)]))




%% saved all step used to produce the final quantification
for io1=1:max(A7(:))
imwrite(uint16(leavs{1,io1}),['Leave_' num2str(io1) '_fmask.tif'],'tif');
imwrite(leavs{2,io1},['Leave_' num2str(io1) '_ori.tif'],'tif');
imwrite(leavs{3,io1},['Leave_' num2str(io1) '_flat.tif'],'tif');
imwrite(uint16(bwlabel(leavs{4,io1})),['Leave_' num2str(io1) '_dmask.tif'],'tif');
% pause
end
