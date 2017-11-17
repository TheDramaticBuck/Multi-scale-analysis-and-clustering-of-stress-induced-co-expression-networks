% AMI.m

% Adjusted Mutual Information function

% Based on original code by Nguyen Xuan Vinh 2011
% (https://uk.mathworks.com/matlabcentral/fileexchange/33144-the-adjusted-mutual-information?requestedDomain=www.mathworks.com)

% Inputs:

% P1 and P2 ...................... clustering solutions for different experiments at the same or at different Markov times (resolution)

% Output:

% AMI_ ........................... Adjusted Mutual Information value for the pair P1:P2

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

function [AMI_]=AMI(P1,P2)

% Build the contingency table from membership arrays
   
R=max(P1);

C=max(P2);

n=length(P2);
 
N=n;

% Identify and remove the missing labels
   
list_t=ismember(1:R,P1);
    
list_m=ismember(1:C,P2);

T=Contingency(P1,P2);

T=T(list_t,list_m);


% Update the true dimensions

[R,C]=size(T);

if C>1 
    a=sum(T');
else
    a=T';
end

if R>1 
    b=sum(T);
else
    b=T;
end

% Calculating the Entropies

Ha=-(a/n)*log(a/n)'; 

Hb=-(b/n)*log(b/n)';

% Calculate the MI (unadjusted for chance)

MI=0;

for i=1:R
    for j=1:C
        if T(i,j)>0 
            MI=MI+T(i,j)*log(T(i,j)*n/(a(i)*b(j)));
        end
    end
end

MI=MI/n;

% Adjust for chance

AB=a'*b;

bound=zeros(R,C);

sumPnij=0;

E3=(AB/n^2).*log(AB/n^2);

EPLNP=zeros(R,C);

LogNij=log((1:min(max(a),max(b)))/N);

for i=1:R
    for j=1:C
        
        sumPnij=0;
        
        nij=max(1,a(i)+b(j)-N);
        
        X=sort([nij N-a(i)-b(j)+nij]);
        
        if N-b(j)>X(2)
            nom=[a(i)-nij+1:a(i) b(j)-nij+1:b(j) X(2)+1:N-b(j)];
            dem=[N-a(i)+1:N 1:X(1)];
        else
            nom=[a(i)-nij+1:a(i) b(j)-nij+1:b(j)];       
            dem=[N-a(i)+1:N N-b(j)+1:X(2) 1:X(1)];
        end
        
        p0=prod(nom./dem)/N;
        
        sumPnij=p0;
        
        EPLNP(i,j)=nij*LogNij(nij)*p0;
        
        p1=p0*(a(i)-nij)*(b(j)-nij)/(nij+1)/(N-a(i)-b(j)+nij+1);  
        
        for nij=max(1,a(i)+b(j)-N)+1:1:min(a(i), b(j))
            
            sumPnij=sumPnij+p1;
            EPLNP(i,j)=EPLNP(i,j)+nij*LogNij(nij)*p1;
            p1=p1*(a(i)-nij)*(b(j)-nij)/(nij+1)/(N-a(i)-b(j)+nij+1);            
            
        end
        
         CC=N*(a(i)-1)*(b(j)-1)/a(i)/b(j)/(N-1)+N/a(i)/b(j);
         
         bound(i,j)=a(i)*b(j)/N^2*log(CC);  
         
    end
end

EMI_bound=sum(sum(bound));

EMI=sum(sum(EPLNP-E3));

AMI_=(MI-EMI)/(max(Ha,Hb)-EMI);

NMI=MI/sqrt(Ha*Hb);


% If expected mutual information is negligible, use NMI.

if abs(EMI)>EMI_bound
    
    fprintf('The EMI is small: EMI < %f, setting AMI=NMI',EMI_bound);
    
    AMI_=NMI;
end

end

    

