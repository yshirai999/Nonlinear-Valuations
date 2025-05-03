function [MDQ] = quantize(MD)
    % Remove outliers
    MD = rmoutliers(MD,'percentiles',[5,95]);
 
    % Parameters
    L = 2^4;
    rng('default'); rng(1);
    
    % Quantize Datasets
    [MDQ, ~, ~] = vqsplit(MD',L);
    MDQ = MDQ';
end