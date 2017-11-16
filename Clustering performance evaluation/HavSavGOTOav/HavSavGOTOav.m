% This program computes the value of linear performance functions between gene expression 
% profiles grouped in clusters
%
% Clustering solutions were computed through stability analysis
%
% Clustering performance functions: Hav (average homogeneity), Sav (average
% separation) and GOTOav (average gene ontology term overlap) at each
% resolution, here quantified by the Markov Time
%
%
% For specific details on the underlying methods, including stability analysis, and the datasets see:
% 
%        https://arxiv.org/abs/1703.02872 
%
%
% Copyright (c) 2017 Nuno R. Nene, University of Cambridge

% Email: nunonene@gmail.com

function [Hav_MT,Sav_MT,GOTOavbp_MT,GOTOavmf_MT,GOTOavcc_MT]=HavSavGOTO(T,probenames,datamatrix,SGDmap_bp,SGDmap_mf,SGDmap_cc,C)

% Inputs:

% T .......................................... Markov Times for each partition solution determined through stability analysis
% probenames ................................. Identifiers for each row of datamatrix
% datamatrix ................................. Row-wise normalized expression matrix: n probenames by p time-points
% SGDmap_bp, SGDmap_mf and SGDmap_cc ......... Container map of Gene Ontology terms for S.cerevisiae (SGD); import with load('SC_CompleteGOterms.mat')
% C .......................................... n probes by length(T) matrix of partition solutions identified by stability analysis; each column is a vector where entries are identified by their respective cluster or community at the respective Markov time

% Outputs:

% Hav_MT ..................................... Average homegeneity for each Markov time in T
% Sav_MT ..................................... Average separation for each Markov time in T
% GOTObp_MT .................................. Average GOTO (bp) for each Markov time in T
% GOTOmf_MT .................................. Average GOTO (mf) for each Markov time in T
% GOTOcc_MT .................................. Average GOTO (cc) for each Markov time in T

% For parfor

poolobj = gcp('nocreate'); % If no pool, do not create new one.

if isempty(poolobj)
    
    poolsize = 0;
    
else
    
    parpool('local',8)
end

% **************************************************************************
% Distance matrix
% **************************************************************************

% Matrix D is the distance between genes according to their tim-dependent expression profile

distance='cosine'; % here we resort to a linear dissimilarity function

D=1-squareform(pdist(datamatrix, distance));

%Assert D is a square matrix of size (ngenes X ngenes) 

assert(size(D,1) == size(datamatrix,1));

% *************************************************************************
% Identify GO terms in each cluster
% Calculate Hav's (Homogeneity (average) based on the chosen metric)
% Calculate Sav's (Separation (average) between clusters)
% Calculate GOTO's(GO Term Overlap)
% *************************************************************************
 
 Hav_MT=zeros(numel(T),1); % for the clustering solution at a specific MT
 Sav_MT=zeros(numel(T),1); 
 
 GOTOavbp_MT=zeros(numel(T),1); % biological process
 GOTOavmf_MT=zeros(numel(T),1); % molecular function
 GOTOavcc_MT=zeros(numel(T),1); % cellular compartment     

 all_probes=1:numel(D(:,1));
 
 for i=1:length(T)
    
%   T_str=num2str(T(i),'%f');
T_str=num2str(i,'%d');
    
  part_unique=unique(C(:,i));%identify genes belonging to each cluster
  
  numPrbsclust=zeros(numel(part_unique),1);
    
  part_order=zeros(numel(part_unique),1);
        
  % sort clusters in ascending Prb number
  
  part_order=zeros(numel(part_unique),1);
     
         
  
  for oj=1:length(part_unique) 
      [r,c] = find(C(:,i)==part_unique(oj)); % for each partition
    
      part_order(oj,1)=numel(r);
  end
  
  
  part_unique=sortrows([part_order part_unique]);
  part_unique=part_unique(:,2);
  
  clear part_order;
         
         
                
  GOTOclust_bp=zeros(length(part_unique),1);% for the GOTO index in each cluster (biological process)
  GOTOclust_mf=zeros(length(part_unique),1);% for the GOTO index in each cluster (molecular function)
  GOTOclust_cc=zeros(length(part_unique),1);% for the GOTO index in each cluster (cellular compartment)
  
  Havclust=zeros(length(part_unique),1); % for the Hav index in each cluster
  Savclust=zeros(length(part_unique),1);% for the Sav index in each cluster
            
    
        
         
        for j=1:length(part_unique)
    
             [r,c] = find(C(:,i)==part_unique(j));
        
         
             numPrbsclust(j,1)=numel(r);% numb probes in each cluster (to be used in the plots)
        
             Hmg=[];
             Sav=[];
             
                % *********************************************************
                % Calculate cluster homogeneity and separation from others
                % *********************************************************
                
                
                if(length(r)>1 && length(r)<length(datamatrix(:,1)))
               
                    Hmg=D(r,r);
                    not_r=setdiff(all_probes,r);
                    Sav=D(r,not_r);
                        
                                     
                    [row,col]=find(triu(Hmg,1));
                    Havclust(j,1)=(sum(sum(triu(Hmg,1),2)))./length(row);
                    Savclust(j,1)=sum(sum(Sav,2))./(numel(Sav(:,1))*numel(Sav(1,:)));
                    
                    
                     
                    % ****************************************************
                    % Find GO terms for each cluster and calculate GOTO
                    % score
                    % ****************************************************
                
                    probes=probenames(r); % see above 
                                 
                 
                    GOterms_withancst_bp=[];
                    geneIndMatrix_withancst_bp=[];
                 
                    GOterms_withancst_mf=[];
                    geneIndMatrix_withancst_mf=[];
                    
                    GOterms_withancst_cc=[];
                    geneIndMatrix_withancst_cc=[];
               
                    for iii = 1:numel(probes)
                        
                    %-----------------------------    
                    %Biological Process GOTO score
                    %-----------------------------
                     if isKey(SGDmap_bp,probes{iii})
                      
                         [geneIndMatrix_withancst_bp,GOterms_withancst_bp]=GOTObp(iii,probes,geneIndMatrix_withancst_bp,GOterms_withancst_bp,SGDmap_bp);
                         
                     else
                      
                          geneIndMatrix_withancst_bp(iii,:)=0;
                     
                     end
                     
                     
                     %-----------------------------
                     %Molecular Function GOTO score
                     %-----------------------------
                     
                     if isKey(SGDmap_mf,probes{iii})
                      
                          [geneIndMatrix_withancst_mf,GOterms_withancst_mf]=GOTOmf(iii,probes,geneIndMatrix_withancst_mf,GOterms_withancst_mf,SGDmap_mf);
                      
                     else
                      
                          geneIndMatrix_withancst_mf(iii,:)=0;
                     
                     end
                     
                    
                    
                    %--------------------------------
                    % Cellular Compartment GOTO score
                    %--------------------------------
                     
                     if isKey(SGDmap_cc,probes{iii})
                      
                       [geneIndMatrix_withancst_cc,GOterms_withancst_cc]=GOTOcc(iii,probes,geneIndMatrix_withancst_cc,GOterms_withancst_cc,SGDmap_cc);
                         
                         
                     else
                      
                          geneIndMatrix_withancst_cc(iii,:)=0;
                     
                     
                     end
                    end
                     
                     
                     
                
                
                
                % *************************************************************
                % Calculate GOTO for clusters at each MT
                % *************************************************************
           
                                
                k=2;
                         
                Cmb = combnk(1:numel(probes),k);% combination 2 and 2 without repetition
            
                
                res_withancst_bp=zeros(numel(Cmb(:,1)),1);%[];
                
                res_withancst_mf=zeros(numel(Cmb(:,1)),1);%[];
                
                res_withancst_cc=zeros(numel(Cmb(:,1)),1);%[];
            
           
            
                
                parfor cmbi=1:numel(Cmb(:,1))
            
                    res_withancst_bp(cmbi,1)=sum((geneIndMatrix_withancst_bp(Cmb(cmbi,1),:)+geneIndMatrix_withancst_bp(Cmb(cmbi,2),:))>1);      
                    res_withancst_mf(cmbi,1)=sum((geneIndMatrix_withancst_mf(Cmb(cmbi,1),:)+geneIndMatrix_withancst_mf(Cmb(cmbi,2),:))>1);      
                    res_withancst_cc(cmbi,1)=sum((geneIndMatrix_withancst_cc(Cmb(cmbi,1),:)+geneIndMatrix_withancst_cc(Cmb(cmbi,2),:))>1);      

                end
                
                res_withancst_bp=sum(res_withancst_bp);
                res_withancst_mf=sum(res_withancst_mf);
                res_withancst_cc=sum(res_withancst_cc);
         
            
                
                GOTOclust_bp(j,1)=(2./(numel(probes)*(numel(probes)-1))).*res_withancst_bp;% GOTO per cluster at a specific MT (bp)
                GOTOclust_mf(j,1)=(2./(numel(probes)*(numel(probes)-1))).*res_withancst_mf;% GOTO per cluster at a specific MT (mf) 
                GOTOclust_cc(j,1)=(2./(numel(probes)*(numel(probes)-1))).*res_withancst_cc;% GOTO per cluster at a specific MT (cc)
            
           
                    
                
                else
                   if(length(r)<length(datamatrix(:,1))) 
                    
                       probes=probenames(r); % see above 
                    
                       Havclust(j,1)=0;
                   
                    
                       not_r=setdiff(all_probes,r);
                    
                       Sav=D(r,not_r);
                    
                       Savclust(j,1)=sum(sum(Sav,2))./(numel(Sav(:,1))*numel(Sav(1,:)));
                    
                    
                       GOTOclust_bp(j,1)=0;
                       GOTOclust_mf(j,1)=0;
                       GOTOclust_cc(j,1)=0;
                   else
                       
                       
                       Hmg=D(r,r);
                       
                       [row,col]=find(triu(Hmg,1));
                    
                       Havclust(j,1)=(sum(sum(triu(Hmg,1),2)))./length(row);
                       
                       Savclust(j,1)=0;
                       
                        
                        
                    % ****************************************************
                    % Find GO terms for each cluster and calculate GOTO
                    % score
                    % ****************************************************
                
                    probes=probenames(r); % see above 
                                 
                 
                    GOterms_withancst_bp=[];
                    geneIndMatrix_withancst_bp=[];
                 
                    GOterms_withancst_mf=[];
                    geneIndMatrix_withancst_mf=[];
                    
                    GOterms_withancst_cc=[];
                    geneIndMatrix_withancst_cc=[];
               
                    for iii = 1:numel(probes)
                        
                    %-----------------------------    
                    %Biological Process GOTO score
                    %-----------------------------
                     if isKey(SGDmap_bp,probes{iii})
                      
                         [geneIndMatrix_withancst_bp,GOterms_withancst_bp]=GOTObp(iii,probes,geneIndMatrix_withancst_bp,GOterms_withancst_bp,SGDmap_bp);
                         
                     else
                      
                          geneIndMatrix_withancst_bp(iii,:)=0;
                     
                     end
                     
                     
                     %-----------------------------
                     %Molecular Function GOTO score
                     %-----------------------------
                     
                     if isKey(SGDmap_mf,probes{iii})
                      
                          [geneIndMatrix_withancst_mf,GOterms_withancst_mf]=GOTOmf(iii,probes,geneIndMatrix_withancst_mf,GOterms_withancst_mf,SGDmap_mf);
                      
                     else
                      
                          geneIndMatrix_withancst_mf(iii,:)=0;
                     
                     end
                     
                    
                    
                    %--------------------------------
                    % Cellular Compartment GOTO score
                    %--------------------------------
                     
                     if isKey(SGDmap_cc,probes{iii})
                      
                       [geneIndMatrix_withancst_cc,GOterms_withancst_cc]=GOTOcc(iii,probes,geneIndMatrix_withancst_cc,GOterms_withancst_cc,SGDmap_cc);
                         
                         
                     else
                      
                          geneIndMatrix_withancst_cc(iii,:)=0;
                     
                     
                     end
                    end
                     
                
                    % *************************************************************
                    % Calculate GOTO for clusters at each MT
                    % *************************************************************
           
                                
                    k=2;
                         
                    Cmb = combnk(1:numel(probes),k);% combination 2 and 2 without repetition
            
                
                    res_withancst_bp=zeros(numel(Cmb(:,1)),1);%[];
                
                    res_withancst_mf=zeros(numel(Cmb(:,1)),1);%[];
                
                    res_withancst_cc=zeros(numel(Cmb(:,1)),1);%[];
            
           
            
                
                    parfor cmbi=1:numel(Cmb(:,1))
           
                        res_withancst_bp(cmbi,1)=sum((geneIndMatrix_withancst_bp(Cmb(cmbi,1),:)+geneIndMatrix_withancst_bp(Cmb(cmbi,2),:))>1);      
                        res_withancst_mf(cmbi,1)=sum((geneIndMatrix_withancst_mf(Cmb(cmbi,1),:)+geneIndMatrix_withancst_mf(Cmb(cmbi,2),:))>1);      
                        res_withancst_cc(cmbi,1)=sum((geneIndMatrix_withancst_cc(Cmb(cmbi,1),:)+geneIndMatrix_withancst_cc(Cmb(cmbi,2),:))>1);      

                    end
                
                    res_withancst_bp=sum(res_withancst_bp);
                    res_withancst_mf=sum(res_withancst_mf);
                    res_withancst_cc=sum(res_withancst_cc);
                    
       
                    GOTOclust_bp(j,1)=(2./(numel(probes)*(numel(probes)-1))).*res_withancst_bp;% GOTO per cluster at a specific MT
                    GOTOclust_mf(j,1)=(2./(numel(probes)*(numel(probes)-1))).*res_withancst_mf;% GOTO per cluster at a specific MT
                    GOTOclust_cc(j,1)=(2./(numel(probes)*(numel(probes)-1))).*res_withancst_cc;% GOTO per cluster at a specific MT
                    
                   end
                      
                end
                
                           
               
                Hav_MT(i,1)=Hav_MT(i,1)+(Havclust(j,1)./numel(part_unique));
                
                Sav_MT(i,1)=Sav_MT(i,1)+(Savclust(j,1)./numel(part_unique));
                  
                GOTOavbp_MT(i,1)=GOTOavbp_MT(i,1)+(numel(probes).*GOTOclust_bp(j,1)./numel(probenames));
                GOTOavmf_MT(i,1)=GOTOavmf_MT(i,1)+(numel(probes).*GOTOclust_mf(j,1)./numel(probenames));
                GOTOavcc_MT(i,1)=GOTOavcc_MT(i,1)+(numel(probes).*GOTOclust_cc(j,1)./numel(probenames));
 
        end
            
 end

    
       

 
    
    
   
     
    
    
    

 
