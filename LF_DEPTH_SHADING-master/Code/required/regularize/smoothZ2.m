function Zsmooth = smoothZ2(Zs, Ws, Gs, lambdas_eq, lambda_flat, lambda_smooth, ROBUSTIFY_SMOOTHNESS, IMAGE, GRADIENT_THRES,SOFTEN_EPSILON,CONVERGE_FRACTION)
% G is the contrast function Gs{1} is for defocus, Gs{2} is for corresp

% gradient analysis 
h           = [0 -1 0; -1 4 -1; 0 -1 0]/8;
im_gradient = imfilter(IMAGE,h,'symmetric');
im_gradient = sqrt((im_gradient(:,:,1).^2+im_gradient(:,:,2).^2+im_gradient(:,:,3).^2)./3);
im_gradient_filt = double(im_gradient < GRADIENT_THRES);

Aeqs = {};
beqs = {};

for Zi = 1:length(Zs)

  Z = Zs{Zi};
  W = Ws{Zi};
  G = Gs{Zi};
  
  Zvalid = ~isnan(Z);
  
  n = numel(Z);
  
  for Gi = 1:size(Gs{1},3)
      
    G_interest = G(:,:,Gi);
    Aeq1 = sparse(1:n, 1:n, W(:));
    Aeq2 = sparse(1:n, 1:n, G_interest(:));
    
    Aeq = Aeq1 .* Aeq2;
    
    Aeq = Aeq(Zvalid,:);
    beq = W(Zvalid) .* reshape(G_interest(Zvalid),n,1) .* Gi;
  
    Aeqs{size(Gs{1},3)*(Zi-1) + Gi} = lambdas_eq(Zi) * Aeq;
    beqs{size(Gs{1},3)*(Zi-1) + Gi} = lambdas_eq(Zi) * beq;
  end
  
end

Aeq = cat(1,Aeqs{:});
beq = cat(1,beqs{:});



if lambda_flat > 0
  f = [1, -1]/2;
  Aflat1 = conv2mat(size(Z), f, im_gradient_filt);
  Aflat2 = conv2mat(size(Z), f', im_gradient_filt);
  Aflat = [Aflat1; Aflat2];
  bflat = sparse(size(Aflat,1), 1);
else
  Aflat = [];
  bflat = [];
end

% if lambda_smooth > 0
%   f = [-1, 2, -1]/4;
%   Asmooth1 = conv2mat(size(Z), f);
%   Asmooth2 = conv2mat(size(Z), f');
%   Asmooth = [Asmooth1; Asmooth2];
%   bsmooth = sparse(size(Asmooth,1), 1);
% else
%   Asmooth = [];
%   bsmooth = [];
% end


% if lambda_flat > 0
%   f = [1, -1; -1, 1];
%   Aflat = conv2mat(size(Z), f);
%   bflat = sparse(size(Aflat,1), 1);
% else
%   Aflat = [];
%   bflat = [];
% end

if lambda_smooth > 0
  f = [0, -1, 0; -1, 4, -1; 0, -1, 0];
  Asmooth = conv2mat(size(Z), f, im_gradient_filt);
  bsmooth = sparse(size(Asmooth,1), 1);
else
  Asmooth = [];
  bsmooth = [];
end


A = [Aeq; lambda_flat * Aflat; lambda_smooth * Asmooth];
b = [beq; lambda_flat * bflat; lambda_smooth * bsmooth];
tic;

X = full(A \ b);

Z_init = reshape(X, size(Zs{1}));
visualizeZ(Z_init); drawnow;

% Z_numer = 0;
% Z_denom = 0;
% 
% for Zi = 1:length(Zs)
%   W = lambdas_eq(Zi) * Ws{Zi};
%   Z_numer = Z_numer + Zs{Zi} .* W;
%   Z_denom = Z_denom + W;
% end
% Z_init = Z_numer ./ Z_denom;
% 
% % Z_init = ones(size(Z_init)) * median(Z_init(:));
% 
% X = Z_init(:);

%SOFTEN_EPSILON = .1;
%SOFTEN_EPSILON = 0.001;
%CONVERGE_FRACTION = 0.001;

losses = [];
for iter = 1:20
  
  delta = abs(A*X-b);
  if ROBUSTIFY_SMOOTHNESS
  
    err = sqrt(delta.^2 + SOFTEN_EPSILON.^2);
    w = 1./err;
    W = sparse(1:length(err), 1:length(err), sqrt(w));

  else
    
    is_robust = false(size(b));
    is_robust(1:length(beq)) = true;
    
    err = zeros(size(delta));
    err(is_robust) = sqrt(delta(is_robust).^2 + SOFTEN_EPSILON.^2);
    err(~is_robust) = 0.5 * delta(~is_robust).^2;
    
    w = 1./err;
    w(~is_robust) = 1;
    W = sparse(1:length(w), 1:length(w), sqrt(w));
    
  end
  loss = mean(err);
  
  fprintf('%02d: %f, %s\n', iter, loss, time2human(toc));
  tic
  
  losses(iter) = loss;
  if (iter > 1) && (abs(losses(iter-1) - losses(iter)) / losses(iter)) < CONVERGE_FRACTION;
    fprintf('Converged\n');
    break
  end
  
  
  WA = W * A;
  Wb = W * b;
  X = full(WA \ Wb);
  
  Zsmooth = reshape(X, size(Z));
  visualizeZ(Zsmooth); drawnow;
  
end


