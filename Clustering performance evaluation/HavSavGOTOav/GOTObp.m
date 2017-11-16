% This program computes the value of the Gene Ontology Term Overlap function
% between genes is the same cluster (biological process)

% For specific details on the underlying methods:
% 
%        https://arxiv.org/abs/1703.02872 
%
%
% Copyright (c) 2017 Nuno R. Nene, University of Cambridge

% Email: nunonene@gmail.com



function [geneIndMatrix_withancst_bp,GOterms_withancst_bp]=GOTObp(iii,probes,geneIndMatrix_withancst_bp,GOterms_withancst_bp,SGDmap_bp)


goid_withancst_bp=SGDmap_bp(probes{iii});
                                                     
goid_withancst_bp=unique(goid_withancst_bp);
                            
                             
if(iii==1)
    
    if(~isempty(goid_withancst_bp))
                         
                                 
        GOterms_withancst_bp=[GOterms_withancst_bp;goid_withancst_bp];
                         
                                 
        geneIndMatrix_withancst_bp(iii,1:length(goid_withancst_bp))=1;
                           
        
    end
    
else
    
    if(~isempty(goid_withancst_bp))
                             
    
        [cia,aia,bia]=intersect(GOterms_withancst_bp,goid_withancst_bp);
        
        [cda,ada]=setdiff(goid_withancst_bp,GOterms_withancst_bp);
                                 
        
        if(~isempty(cia))
        
            geneIndMatrix_withancst_bp(iii,aia)=1;
                                
        else
            
            geneIndMatrix_withancst_bp(iii,:)=0;
                                
        end
        
        
        if(~isempty(cda))
                                    
            geneIndMatrix_withancst_bp=[geneIndMatrix_withancst_bp zeros(iii,length(ada))];
                                     
            geneIndMatrix_withancst_bp(iii,length(GOterms_withancst_bp)+1:end)=1;
                                     
            GOterms_withancst_bp=[GOterms_withancst_bp cda];
                                 
        end
        
        
    end
    
end
                         
                

end