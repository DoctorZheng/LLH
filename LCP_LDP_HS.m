function result = LCP_LDP_HS(image)
    if ndims(image) == 3
        image = rgb2gray(image);
    end
    image = double(image);
    [rows, cols] = size(image);

    if rows < 5 || cols < 5
        error('Input image must be at least 5x5 pixels.');
    end

    valid_r = rows - 4;
    valid_c = cols - 4;

    LCP_map = zeros(valid_r, valid_c);
    LDP_map_D = zeros(valid_r, valid_c); % binary count (0～8)
    LDP_map = zeros(valid_r, valid_c);   % weighted value
    HS_map = zeros(valid_r, valid_c);

    [Ix, Iy] = gradient(image);
    [Ixx, ~] = gradient(Ix);
    [~, Iyy] = gradient(Iy);
    divV = Ixx + Iyy;

    Ic = image(3:end-2, 3:end-2);  % (valid_r x valid_c)

    dirs = [-2,-2; -2,0; -2,2; ...
             0, 2;  2,2;  2,0; ...
             2,-2;  0,-2];  % 8 directions

    N = valid_r * valid_c;
    diffs_LCP = zeros(8, N);
    diffs_LDP = zeros(8, N);

    for k = 1:8
        dx = dirs(k,1); dy = dirs(k,2);
        far = image(3+dx : end-2+dx, 3+dy : end-2+dy);
        near = image(3+dx/2 : end-2+dx/2, 3+dy/2 : end-2+dy/2);
        
        diffs_LCP(k, :) = Ic(:)' + far(:)' - 2*near(:)';
        diffs_LDP(k, :) = (far(:)' - Ic(:)') / 2;
    end

    LCP_vals = sum(diffs_LCP > 0, 1);  % 1 x N
    LCP_map = reshape(LCP_vals, valid_r, valid_c);

    total_binary = sum(diffs_LDP > 0, 1);
    LDP_map_D = reshape(total_binary, valid_r, valid_c);
    sum_diffs_vec = sum(diffs_LDP, 1);  % 1×N

    kernel_8 = ones(3); 
    kernel_8(2,2) = 0;  
    neighbor_sum_full = imfilter(image, kernel_8, 'replicate'); % same size as image
    neighbor_sum_center = neighbor_sum_full(3:end-2, 3:end-2);   % [valid_r × valid_c]
    Iref_map = neighbor_sum_center / 8;               

    Ic_vec = Ic(:);                  % N×1
    Iref_vec = Iref_map(:);          % N×1
    sum_diffs_vec = sum_diffs_vec';  % 转为 N×1

    LDP_vec = zeros(N, 1);

    idx1 = (Ic_vec ~= 0);
    LDP_vec(idx1) = (sum_diffs_vec(idx1) ./ Ic_vec(idx1)) + 4;

    idx2 = (~idx1) & (Iref_vec ~= 0);
    LDP_vec(idx2) = sum_diffs_vec(idx2) ./ Iref_vec(idx2); 

    LDP_map = reshape(LDP_vec, valid_r, valid_c);

    neigh = [-1,-1; -1,0; -1,1; 0,1; 1,1; 1,0; 1,-1; 0,-1];
    HS_vals = zeros(1, N);
    cnt = 1;
    for x = 3:rows-2
        for y = 3:cols-2
            hs_count = 0;
            for k = 1:8
                nx = x + neigh(k,1);
                ny = y + neigh(k,2);
                if divV(nx, ny) < 0
                    hs_count = hs_count + 1;
                end
            end
            HS_vals(cnt) = hs_count;
            cnt = cnt + 1;
        end
    end
    HS_map = reshape(HS_vals, valid_r, valid_c);

    bins_edges = 0:8;  % 9 bins: 0,1,...,8

    h1 = histcounts(LCP_map(:), [bins_edges, Inf]);
    h2 = accumarray(LDP_map_D(:) + 1, LDP_map(:), [9, 1]);  % LDP_map_D ∈ [0,8]
    h3 = histcounts(HS_map(:), [bins_edges, Inf]);

    h1 = h1 / sum(h1 + eps);
    h2 = h2 / sum(h2 + eps);
    h3 = h3 / sum(h3 + eps);

    result = [h1(:)', h2(:)', h3(:)'];  
end
