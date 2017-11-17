% RunAMI.m

% Program for calculating the Adjusted Mutual Information (AMI) between
% partitions calculated by stability analysis

% Inputs:

% C1,C2 ................... Partition solution found by stability analysis at several resolutions (columns); rows correspond to probes or genes 

% Output:

% AMImat .................. Adjusted Mutual Information matrix between all the partitions


% Based on original code by Nguyen Xuan Vinh 2011
% (https://uk.mathworks.com/matlabcentral/fileexchange/33144-the-adjusted-mutual-information?requestedDomain=www.mathworks.com)

% References:

% [1] 'A Novel Approach for Automatic Number of Clusters Detection based on Consensus Clustering', 
%       N.X. Vinh, and Epps, J., in Procs. IEEE Int. Conf. on 
%       Bioinformatics and Bioengineering (Taipei, Taiwan), 2009.
% [2] 'Information Theoretic Measures for Clusterings Comparison: Is a
%	    Correction for Chance Necessary?', N.X. Vinh, Epps, J. and Bailey, J.,
%	    in Procs. the 26th International Conference on Machine Learning (ICML'09)
% [3] 'Information Theoretic Measures for Clusterings Comparison: Variants, Properties, 
%       Normalization and Correction for Chance', N.X. Vinh, Epps, J. and
%       Bailey, J., Journal of Machine Learning Research, 11(Oct), pages
%       2837-2854, 2010



% Nuno R. Nene, University of Cambridge

% Email: nunonene@gmail.com

function [AMImat]=RunAMI(C1,C2)
 
% For parfor

poolobj = gcp('nocreate'); % If no pool, do not create new one.

if isempty(poolobj)
    parpool('local',2) % choose number of workers according to your specifications
end

% Create AMI mat

AMImat=zeros(numel(C1(1,:)),numel(C2(1,:)));

for i=1:numel(C1(1,:))
    
    C=C1(:,i)+1; % AMImat should be symmetric
    
   parfor j=1:numel(C2(1,:))
       
      [AMImat(i,j)]= AMI(C',C2(:,j)'+1);
       
       
   end
   
end

delete(gcp('nocreate'));

end

