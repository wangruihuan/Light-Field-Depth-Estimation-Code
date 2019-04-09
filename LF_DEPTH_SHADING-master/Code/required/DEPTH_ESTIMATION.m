function [depth] = DEPTH_ESTIMATION(response,minmax)
%DEPTH_ESTIMATION
%           Takes the response and finds the min or max
%           Input : response
%                   minmax        (0 for min- corresp, 1 for max - defocus)
%           Output: depth

%           EQUATION (6) in paper


if(minmax == 1)
    [garbage depth] = max(response,[],3)  ;
else
    [garbage depth] = min(response,[],3)  ;
end

% depth = zeros(size(response,1),size(response,2));
% 
% for i = 1:size(response,2)
%     for j = 1:size(response,1)
%         [peaks, locs] = findpeaks(-squeeze(response(j,i,:)),'SortStr','descend');
%         if(numel(locs)~=0)
%             depth(j,i) = locs(1);
%         else
%             depth(j,i) = find(-squeeze(response(j,i,:))==max(-squeeze(response(j,i,:))));
%         end
%     end
% end

end

