%% BLOCK BOOTSTRAP %%
% This function performs block bootstrapping of data
% Reshuffling is done across columns

% INPUTS %
% x = data to be reshuffled
% blocks = block lenght

% OUTPUTS %
% x_b = reshuffled series

function x_b = block(x,blocks); 
for i=1:blocks:(size(x,2)-blocks+1);
    j=randi(size(x,2)-blocks+1);
    x_b(:,i:i+blocks-1) = x(:,j:j+blocks-1); 
end