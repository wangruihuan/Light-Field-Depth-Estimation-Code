x=(1:64)';
cost1=zeros(64,1);
for i=1:1
    for j=1:1
      cost1=squeeze(cost_h(i,j,1:64));
      plot(x,cost1,'-o'); 
      clear cost1;
      hold on;
    end
end
hold off;
nD=64;
save_im=(256/(nD))*(labels_max-1);
