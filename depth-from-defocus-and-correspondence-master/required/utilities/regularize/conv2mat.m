function A = conv2mat(sz, Fm, grad_filt)

F = reshape(Fm(end:-1:1), size(Fm));

[i0, j0] = find(true(sz).*grad_filt);
idx0 = sub2ind(sz, i0, j0);

idxs = {};
fs = {};
for oi = 1:size(F,1)
  for oj = 1:size(F,2)
    
    f = F(oi,oj);
    if f == 0
      continue;
    end
    
    i = i0 + (oi-1);
    j = j0 + (oj-1);
    
    keep = (i <= sz(1)) & (j <= sz(2));
    idx = sub2ind(sz, i(keep), j(keep));
    idxs{end+1} = nan(size(idx0));
    idxs{end}(keep) = idx;
    fs{end+1} = f;
%     m = length(idx);
%     sparse(1:m, idx0, 1, m, n) - sparse(1:m, idx_curve(:,1), 1, m, n) - sparse(1:m, idx_curve(:,3), 1, m, n)
%     [idx0(keep), idx, F(oi,oj)
  end
end
idxs = cat(2,idxs{:});
idxs = idxs(all(~isnan(idxs),2),:);

n = prod(sz);
m = size(idxs,1);
A = sparse(0);
for j = 1:size(idxs,2)
  A = A + sparse(1:m, idxs(:,j), fs{j}, m, n);
end


end