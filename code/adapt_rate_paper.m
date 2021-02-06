


clear;
clc;

% initialization
M = 400; %train images number
NUM = 200; %test images number
MAXN = 1000;


r = 0.3; % upper bound of sampling rate
G_var = 0.01;%noise varience


% random vector for selecting train data and test data
R = load('R.mat');
R = R.R; 

N = int32(256*256*(r));

%upper bond of performance perservation
u_a = 0.8; 


m=1:256;
n=1:256;
z=power(-1,(bsxfun(@plus,m,n')));
% real and image
a = real(dftmtx(256));
b = imag(dftmtx(256));
p_o = zeros(256,256);
p_n = zeros(256,256);

%load data
f_data = 0; % 1= ibsr
if f_data ==1
        fpath = 't1_IBSR';
        flist = dir(sprintf('%s/*.png', fpath));
        images = [];
        for imidx = 1:min(length(flist), 1000)
            fprintf('[%d]', imidx); 
            fname = sprintf('%s/%s', fpath, flist(imidx).name);
            im = imread(fname);
            im = rgb2gray(im);
            images{length(images)+1} = im;
        end
        fprintf('\n');
else
        fpath = 't1_MGC'; 
        flist = dir(sprintf('%s/*.png', fpath));

        images = [];
        for imidx = 1:min(length(flist), 1000)
            fprintf('[%d]', imidx);
            fname = sprintf('%s/%s', fpath, flist(imidx).name);
            im = imread(fname);
            %im = rgb2gray(im);
            images{length(images)+1} = im;
        end
        %imagesc(im);
        fprintf('\n');
end

%generate noise data
for j = 1:MAXN
        P = images{R(j)};
        P=im2double(P);
        noise_r = G_var.*randn(256);
        noise_i = G_var.*randn(256);
        P = (P+noise_r).*(P+noise_r) + noise_i.*noise_i;
        P = sqrt(P);
        noise_images{j} = P;
end








% train 
for j = 1:M
    P = noise_images{j};
    P0 = images{R(j)};
    P=im2double(P);
    P0=im2double(P0);
    P = P.* z;
    P0 = P0.* z;
    c = a'*P*a - b'*P*b;
    d = a'*P0*a - b'*P0*b;
    e = a'*P *b + b'*P *a;
    f = a'*P0 *b + b'*P0 *a;
    p_n = p_n + (c-d).^2 + (e-f).^2;
   % p_n =  p_n + c.^2  + e.^2 -2 * (c.*d + e.*f);
    p_o = p_o + f.^2 + d.^2;
end
p_k = p_n - p_o;


r_num = 0;
[B2,index2] = sort(p_k(:),'ascend');
p3 = zeros(256,256);
[rc,c]=ind2sub(size(p_k),index2);
for k = 1: 256*256
    if  p_k(rc(k),c(k)) <= 0
    r_num = r_num +1;
    end 
end
r_s = r_num/(256*256)


if r <= r_s
    flag = 1;
    for k = 1:N
        p3(rc(k),c(k)) = 1;
    end
else 
    flag = 0;
 for k = 1:int32(256*256*(r_s))
        p3(rc(k),c(k)) = 1;
    end
end

%sampling regualization

n = 1:M;
for m = 1:M
    img = noise_images{m};
    Y = fft2(img);
    Y = fftshift(Y);
    y2 = Y.*p3;
   % y = ifft2(y2);
    P=im2double(images{R(m)});
    P = fft2(P);
    P = fftshift(P);
    n(m) =  sum(sum(abs((P-y2).^2)))/(256*256);
end
n_0 = mean(n);

if u_a > 0 
    low = 1; 
    if r_s > r
        high = int32(256*256*(r));
    else
        high = int32(256*256*(r_s));
    end
    flag_num = 0;
    flag_num_b = 0;
    while low <= high 
        flag_num = flag_num +1;
         mid = fix((low + high)/2);
         p3 = zeros(256,256);
        for k = 1:mid
            p3(rc(k),c(k)) = 1;
        end
        n = 1:M;
        for m = 1:M
            img = noise_images{m};
            Y = fft2(img);
            Y = fftshift(Y);
            y2 = Y.*p3;
           % y = ifft2(y2);
            P=im2double(images{R(m)});
            P = fft2(P);
            P = fftshift(P);
            n(m) =  sum(sum(abs((P-y2).^2)))/(256*256);
        end
        n_1 = mean(n);
        if n_1 - n_0 > u_a
            p3 = zeros(256,256);
            for k = 1:mid+1
                p3(rc(k),c(k)) = 1;
            end
            n = 1:M;
            for m = 1:M
                img = noise_images{m};
                Y = fft2(img);
                Y = fftshift(Y);
                y2 = Y.*p3;
                P=im2double(images{R(m)});
                P = fft2(P);
                P = fftshift(P);
                n(m) =  sum(sum(abs((P-y2).^2)))/(256*256);
            end
            n_2 = mean(n);
        else
            p3 = zeros(256,256);
            for k = 1:mid-1
            p3(rc(k),c(k)) = 1;
            end
            n = 1:M;
            for m = 1:M
                img = noise_images{m};
                Y = fft2(img);
                Y = fftshift(Y);
                y2 = Y.*p3;
                P=im2double(images{R(m)});
                P = fft2(P);
                P = fftshift(P);
                n(m) = sum(sum(abs((P-y2).^2)))/(256*256);
            end
            n_2 = mean(n);
        end
        
        if (n_2-n_0- u_a)*(n_1-n_0 - u_a) < 0 
            result_index = mid;
            flag_num_b = 1;
            break ;                        
        end
        if (n_1 - n_0 -u_a) > 0  
             low = mid + 1;
             n_cnum = n_1 - n_0 -u_a;
             n_cnum_b = (n_2-n_0- u_a)*(n_1-n_0 - u_a);
        else                      
             high = mid - 1;
             n_enum = n_1 - n_0 -u_a;
             n_enum_b = (n_2-n_0- u_a)*(n_1-n_0 - u_a);
        end
        if low > high              
            result_index = NaN;
            flag_num_b = 2;
            break ;
        end
    end
    learned_rate = double(mid)/(256*256)%learned sampling rate

    p3 = zeros(256,256);
    for k = 1:int32(256*256*(learned_rate))
        p3(rc(k),c(k)) = 1; % measurement matrix 
    end
end




% test
n = M+1:M+1+NUM;
ss = M+1:M+1+NUM;
NM_n = M+1:M+1+NUM;
rlne = M+1:M+1+NUM;
for m = M+1:M+1+NUM
    img = noise_images{m};
    Y = fft2(img);
    Y = fftshift(Y);
    y2 = Y.*p3;
    y = ifft2(y2);
    P=im2double(images{R(m)});
    n(m-M) = psnr(P,abs(y));
    ss(m-M) = ssim(P,abs(y));
    NM_n(m-M) = nmse_fun(P,abs(y));
end
% n
v_psnr = mean(n)%PSNR
v_ssim = mean(ss)%SSIM
v_nmse = mean(NM_n)*100%NMSE










