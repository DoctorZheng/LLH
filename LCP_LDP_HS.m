function result = LCP_LDP_HS(image)
    
[rows, cols, h] = size(image); % 获取图像的尺寸
if h==3
    image =  rgb2gray(image);
end
image=double(image);

LCP_map = zeros(rows-4, cols-4); % 初始化LDP结果图像
LDP_map_D = zeros(rows-4, cols-4); % 初始化LDP结果图像
LDP_map = zeros(rows-4, cols-4); % 初始化LDP结果图像
HS_map = zeros(rows-4, cols-4); % 初始化LDP结果图像

% 计算图像的梯度
[Ix, Iy] = gradient(image); 
% 计算梯度的散度（哈密顿散度）
divV = Ix + Iy;

% 遍历每个像素点（跳过边界像素）
for x = 3:rows-2
    for y = 3:cols-2
        % 获取中心像素 I_c
        I_c = image(x, y);
        % 获取5×5窗口中的中心像素的8个相邻像素点的坐标
        neighbors = [
            x-1, y-1; x-1, y; x-1, y+1;
            x, y+1; x+1, y+1; x+1, y; x+1, y-1;
            x, y-1
        ]; 

        %定义星状射线的8个方向 向8个方向计算灰度中心差分?I_j
        diffs = [
            I_c + image(x-2, y-2) - 2*image(x-1, y-1); % 上左
            I_c + image(x-2, y) - 2*image(x-1, y);   % 上
            I_c + image(x-2, y+2) - 2*image(x-1, y+1); % 上右
            I_c + image(x, y+2) - 2*image(x, y+1);   % 右
            I_c + image(x+2, y+2) - 2*image(x+1, y+1); % 下右
            I_c + image(x+2, y) - 2*image(x+1, y);   % 下
            I_c + image(x+2, y-2) - 2*image(x+1, y-1); % 下左
            I_c + image(x, y-2) -  2*image(x, y-1)  % 左
        ];
        con = diffs > 0;        
        % 计算总的灰度变化
        LCP = sum(con);          
        %存储当前像素的LDP值
        LCP_map(x-2, y-2) = LCP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % 定义星状射线的8个方向 向8个方向计算灰度中心差分
        diffs = [
            image(x-2, y-2) - I_c; % 上左
            image(x-2, y) - I_c;   % 上
            image(x-2, y+2) - I_c; % 上右
            image(x, y+2) - I_c;   % 右
            image(x+2, y+2) - I_c; % 下右
            image(x+2, y) - I_c;   % 下
            image(x+2, y-2) - I_c; % 下左
            image(x, y-2) - I_c  % 左
        ]/2;
        total_diff = sum(diffs>0);     
        if I_c ~= 0  %-4~4*255
            LDP = sum(diffs)/I_c + 4;%为了保证后续归一化，将值+4保证正数。
        else
            LDP = 256;
        end 
        % 存储当前像素的LDP值
        LDP_map_D(x-2, y-2) = total_diff;
        LDP_map(x-2, y-2) = LDP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5        
        % 初始化哈密顿散度的聚合（发散）数量
        HS = 0;
        for i = 1:size(neighbors, 1)           
            if divV(neighbors(i,1),neighbors(i,2)) < 0
                HS = HS + 1; % 聚合
            end
        end
        % 将HS存储在HS_map中
        HS_map(x-2, y-2) = HS;         
     end
end    

bins= 8; 
result(1,:)=hist(LCP_map(:),0:bins);
result(2,:) = accumarray(LDP_map_D(:)+1, LDP_map(:), [9 1]);
result(3,:)=hist(HS_map(:),0:bins);
% %归一化到0-1：
result(1,:)=result(1,:)/sum(result(1,:));
result(2,:)=result(2,:)/sum(result(2,:));
result(3,:)=result(3,:)/sum(result(3,:));