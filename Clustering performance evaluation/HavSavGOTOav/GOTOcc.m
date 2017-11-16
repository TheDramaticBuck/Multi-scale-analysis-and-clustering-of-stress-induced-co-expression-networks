% This program computes the value of the Gene Ontology Term Overlap function
% between genes is the same cluster (cellular compartment)

% For specific details on the underlying methods:
% 
%        https://arxiv.org/abs/1703.02872 
%
%
% Copyright (c) 2017 Nuno R. Nene, University of Cambridge

% Email: nunonene@gmail.com

function [geneIndMatrix_withancst_cc,GOterms_withancst_cc]=GOTOcc(iii,probes,geneIndMatrix_withancst_cc,GOterms_withancst_cc,SGDmap_cc)


goid_withancst_cc=SGDmap_cc(probes{iii});
                                                     
goid_withancst_cc=unique(goid_withancst_cc);
                            
                             
if(iii==1)
    
    if(~isempty(goid_withancst_cc))
                         
                                 
        GOterms_withancst_cc=[GOterms_withancst_cc;goid_withancst_cc];
                         
                                 
        geneIndMatrix_withancst_cc(iii,1:length(goid_withancst_cc))=1;
                           
        
    end
    
else
    
    if(~isempty(goid_withancst_cc))
                             
    
        [cia,aia,bia]=intersect(GOterms_withancst_cc,goid_withancst_cc);
        
        [cda,ada]=setdiff(goid_withancst_cc,GOterms_withancst_cc);
                                 
        
        if(~isempty(cia))
        
            geneIndMatrix_withancst_cc(iii,aia)=1;
                                
        else
            
            geneIndMatrix_withancst_cc(iii,:)=0;
                                
        end
        
        
        if(~isempty(cda))
                                    
            geneIndMatrix_withancst_cc=[geneIndMatrix_withancst_cc zeros(iii,length(ada))];
                                     
            geneIndMatrix_withancst_cc(iii,length(GOterms_withancst_cc)+1:end)=1;
                                     
            GOterms_withancst_cc=[GOterms_withancst_cc cda];
                                 
        end
        
        
    end
    
end
                         
                

end