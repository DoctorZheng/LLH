%针对10*10区域的所有纹理的主函数
close all;
clear all;


file_path =  '.\1\';% 图像文件夹路径，
cc=colormap(lines(100));
img_path_list = dir(strcat(file_path,'*.bmp'));%获取该文件夹中所有jpg格式的图像。
img_num = length(img_path_list);%获取图像总数量，
NUMD=zeros(img_num,10);
Dis1=[];
step=4;
for ii = 1:10 %逐一读取图像，    
    I=imread(strcat(num2str(ii),'.bmp'));
        image_size=size(I);
        dimension=numel(image_size);
        if dimension~=2
           IT1 = rgb2gray(I);
        else
           IT1 = (I);             
        end
    [H W]=size(IT1);    
    if ii == 1
        Dis1 = zeros(floor(H*W/(step*step)),img_num, 10);  
    end
    
    IT1 = imgaussfilt(IT1, 2); 
    
    file_path =  strcat(strcat('.\',num2str(ii)),'\');
    if img_num > 0 
        Dis=zeros(floor(H*W/(step*step)),img_num);
        for jj = 1:img_num %逐一读取图像，
            image_name = img_path_list(jj).name;% 图像名，
            I1 =  imread(strcat(file_path,image_name));
            
            image_size=size(I1);
            dimension=numel(image_size);
            if dimension~=2
               IT2 = rgb2gray(I1);
            else
               IT2 = (I1);                
            end
            IT2 = imgaussfilt(IT2, 2); 

            Dnum=0;  
            nuuu=1;  
            
            sized=3;
            for i = sized+1 : step : H-sized
                for j = sized+1 : step : W-sized                                    
                    IT1=I(i-sized:i+sized, j-sized:j+sized);
                    IT2=I1(i-sized:i+sized, j-sized:j+sized);

                    H1= LCP_LDP_HS(IT1);                    
                    H2= LCP_LDP_HS(IT2);                      
                    Dis(nuuu,jj) = pdist2(H1(:)',H2(:)','euclidean'); 
                    if  Dis(nuuu,jj)>1  
                       Dnum=Dnum+1;
                    end  
                    
                     nuuu=nuuu+1;     
                end
            end            

            NUMD(jj,ii)=Dnum/(nuuu);
        end
    end


TT=0.8;
AA=NUMD;
AA(NUMD<TT)=1;
AA(NUMD>=TT)=0;
AA

