close all
image = imread(['EXP6_Arabidopsis_C2_aligment_221111.png']); % fluoscence image that has been alinged to histology image
metadata=importdata(['EXP6_Arabidopsis_C2_alignment.tsv']); % output file of visium aligning program that provides locations for each st spot
%%
imageSize = size(image);
position = metadata.data(:,[6 5]);
coordinate=metadata.data(:,1:2);
num_spot=length(metadata.data(:,1));

[~,dist]=knnsearch(position,position,'k',2);
r=round(mean(dist(:,2))*(55/100/2));
[x,y]=meshgrid((-r:r)/r,(-r:r)/r);
se=abs(x+y*1i)<=1;

%% compute isTissue
indexSpot=find(boolean(metadata.data(:,7)));
position = position(indexSpot,:);
coordinate = coordinate(indexSpot,:);

position = round(position);
isTissue=sparse(position(:,1),position(:,2),ones(size(position(:,1))),imageSize(1),imageSize(2));
isTissue=full(boolean(isTissue));
isTissue=imdilate(isTissue,se);
indexTissue=find(isTissue);

%% compute identitySpot
[Y,X]=meshgrid(1:imageSize(2),1:imageSize(1));
identitySpot=knnsearch(position,[X(isTissue) Y(isTissue)]);

%% compute IndexNeighbor
data=cell(size(position,1),1);
image=rgb2gray(image);
image=double(image);
bw=image~=255;
bw=imfill(bw,"holes");
image(~bw)=0;
value=zeros(num_spot,1);
image=zeros(size(image));
for j=1:size(position,1)
    data{j}=indexTissue(identitySpot==j);
    value(indexSpot(j),1)=max(image(data{j}));
    value(indexSpot(j),2)=mean(image(data{j}));
    image(data{j})=value(indexSpot(j),1);
end
%% save result

f=fopen('EXP6_Arabidopsis_C2_aligment_221111_new.csv','w');
fprintf(f,'x y max mean\n');
for j=1:size(coordinate,1)
    fprintf(f,'%d %d %d %.4f\n',coordinate(j,1),coordinate(j,2),value(j,1),value(j,2));
end
fclose(f);

barcode=barcode(orderAlphabet,:);
data=data(orderAlphabet,:);
figure
imshow(image,[]);
imwrite(uint8(image),'EXP6_Arabidopsis_C2_aligment_221111_max.png')
