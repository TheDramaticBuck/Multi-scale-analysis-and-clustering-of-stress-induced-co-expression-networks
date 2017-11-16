% This program computes the value of the Gene Ontology Term Overlap function
% between genes is the same cluster (molecular function)

% For specific details on the underlying methods:
% 
%        https://arxiv.org/abs/1703.02872 
%
%
% Copyright (c) 2017 Nuno R. Nene, University of Cambridge

% Email: nunonene@gmail.com

function [geneIndMatrix_withancst_mf,GOterms_withancst_mf]=GOTOcc(iii,probes,geneIndMatrix_withancst_mf,GOterms_withancst_mf,SGDmap_mf)


goid_withancst_mf=SGDmap_mf(probes{iii});
                                                     
goid_withancst_mf=unique(goid_withancst_mf);
                            
                             
if(iii==1)
    
    if(~isempty(goid_withancst_mf))
                         
                                 
        GOterms_withancst_mf=[GOterms_withancst_mf;goid_withancst_mf];
                         
                                 
        geneIndMatrix_withancst_mf(iii,1:length(goid_withancst_mf))=1;
                           
        
    end
    
else
    
    if(~isempty(goid_withancst_mf))
                             
    
        [cia,aia,bia]=intersect(GOterms_withancst_mf,goid_withancst_mf);
        
        [cda,ada]=setdiff(goid_withancst_mf,GOterms_withancst_mf);
                                 
        
        if(~isempty(cia))
        
            geneIndMatrix_withancst_mf(iii,aia)=1;
                                
        else
            
            geneIndMatrix_withancst_mf(iii,:)=0;
                                
        end
        
        
        if(~isempty(cda))
                                    
            geneIndMatrix_withancst_mf=[geneIndMatrix_withancst_mf zeros(iii,length(ada))];
                                     
            geneIndMatrix_withancst_mf(iii,length(GOterms_withancst_mf)+1:end)=1;
                                     
            GOterms_withancst_mf=[GOterms_withancst_mf cda];
                                 
        end
        
        
    end
    
end
                         
                

end