function result = LCP_LDP_HS(image)
    
[rows, cols, h] = size(image); % ��ȡͼ��ĳߴ�
if h==3
    image =  rgb2gray(image);
end
image=double(image);

LCP_map = zeros(rows-4, cols-4); % ��ʼ��LDP���ͼ��
LDP_map_D = zeros(rows-4, cols-4); % ��ʼ��LDP���ͼ��
LDP_map = zeros(rows-4, cols-4); % ��ʼ��LDP���ͼ��
HS_map = zeros(rows-4, cols-4); % ��ʼ��LDP���ͼ��

% ����ͼ����ݶ�
[Ix, Iy] = gradient(image); 
% �����ݶȵ�ɢ�ȣ����ܶ�ɢ�ȣ�
divV = Ix + Iy;

% ����ÿ�����ص㣨�����߽����أ�
for x = 3:rows-2
    for y = 3:cols-2
        % ��ȡ�������� I_c
        I_c = image(x, y);
        % ��ȡ5��5�����е��������ص�8���������ص������
        neighbors = [
            x-1, y-1; x-1, y; x-1, y+1;
            x, y+1; x+1, y+1; x+1, y; x+1, y-1;
            x, y-1
        ]; 

        %������״���ߵ�8������ ��8���������Ҷ����Ĳ��?I_j
        diffs = [
            I_c + image(x-2, y-2) - 2*image(x-1, y-1); % ����
            I_c + image(x-2, y) - 2*image(x-1, y);   % ��
            I_c + image(x-2, y+2) - 2*image(x-1, y+1); % ����
            I_c + image(x, y+2) - 2*image(x, y+1);   % ��
            I_c + image(x+2, y+2) - 2*image(x+1, y+1); % ����
            I_c + image(x+2, y) - 2*image(x+1, y);   % ��
            I_c + image(x+2, y-2) - 2*image(x+1, y-1); % ����
            I_c + image(x, y-2) -  2*image(x, y-1)  % ��
        ];
        con = diffs > 0;        
        % �����ܵĻҶȱ仯
        LCP = sum(con);          
        %�洢��ǰ���ص�LDPֵ
        LCP_map(x-2, y-2) = LCP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        % ������״���ߵ�8������ ��8���������Ҷ����Ĳ��
        diffs = [
            image(x-2, y-2) - I_c; % ����
            image(x-2, y) - I_c;   % ��
            image(x-2, y+2) - I_c; % ����
            image(x, y+2) - I_c;   % ��
            image(x+2, y+2) - I_c; % ����
            image(x+2, y) - I_c;   % ��
            image(x+2, y-2) - I_c; % ����
            image(x, y-2) - I_c  % ��
        ]/2;
        total_diff = sum(diffs>0);     
        if I_c ~= 0  %-4~4*255
            LDP = sum(diffs)/I_c + 4;%Ϊ�˱�֤������һ������ֵ+4��֤������
        else
            LDP = 256;
        end 
        % �洢��ǰ���ص�LDPֵ
        LDP_map_D(x-2, y-2) = total_diff;
        LDP_map(x-2, y-2) = LDP;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5        
        % ��ʼ�����ܶ�ɢ�ȵľۺϣ���ɢ������
        HS = 0;
        for i = 1:size(neighbors, 1)           
            if divV(neighbors(i,1),neighbors(i,2)) < 0
                HS = HS + 1; % �ۺ�
            end
        end
        % ��HS�洢��HS_map��
        HS_map(x-2, y-2) = HS;         
     end
end    

bins= 8; 
result(1,:)=hist(LCP_map(:),0:bins);
result(2,:) = accumarray(LDP_map_D(:)+1, LDP_map(:), [9 1]);
result(3,:)=hist(HS_map(:),0:bins);
% %��һ����0-1��
result(1,:)=result(1,:)/sum(result(1,:));
result(2,:)=result(2,:)/sum(result(2,:));
result(3,:)=result(3,:)/sum(result(3,:));