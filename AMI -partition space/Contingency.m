% Contigency function

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


function Cont=Contingency(Mem1,Mem2)

if nargin < 2
   
    error('Contingency: Requires two vector arguments')
  
end

Cont=zeros(max(Mem1),max(Mem2));

for i = 1:length(Mem1)
   Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
end

end       