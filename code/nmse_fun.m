% function nmse = nmse_fun(orgSig,recSig,varargin)
% if isempty(varargin)
%     boun = 0;
% else boun = varargin{1};
% end
% if size(orgSig,2)==1       % if signal is 1-D
%     orgSig = orgSig(boun+1:end-boun,:);
%     recSig = recSig(boun+1:end-boun,:);
% else                       % if signal is 2-D or 3-D
%     orgSig = orgSig(boun+1:end-boun,boun+1:end-boun,:);
%     recSig = recSig(boun+1:end-boun,boun+1:end-boun,:);
% end
% mse=norm(orgSig(:)-recSig(:),2)^2/length(orgSig(:));
% sigEner=norm(orgSig(:))^2;
% nmse=(mse/sigEner);


function nmse = nmse_fun(org_pic,cop_pic)
len = length(org_pic);

fenzi = 0;
fenmu = 0;
for i = 1:1:len
    for j = 1:1:len
        fenzi = fenzi + (org_pic(i,j) - cop_pic(i,j))^2;
        fenmu = fenmu + (org_pic(i,j)^2);
    end
end
nmse = fenzi / fenmu ;