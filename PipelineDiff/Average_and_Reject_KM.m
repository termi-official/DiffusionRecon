function [Dcm2 enum2]= Average_and_Reject_KM(Dcm, enum,ADC_limit)

% Perform averaging on the DWI and nDWI matrices
% 
% Update the number of averages in enum
%
%
% SYNTAX:  [Dcm2 enum2]= Average_and_Reject_KM(Dcm, enum);
%  
%
% INPUTS:   Dcm - DWI image matrix
%                 [y x slices b-values directions averages dataset]
%           
%           enum - Structure which contains information about the dataset 
%          
% OUTPUTS:  Dcm2 - DWI image matrix 
%                 [y x slices b-values directions averages dataset]
%           
%           enum2 - Structure which contains information about the dataset 
%
% Kevin Moulin 08.31.2017
% Kevin.Moulin.26@gmail.com
% Ennis Lab @ UCLA; http://mrrl.ucla.edu
 
    enum2=enum; 
    Dcm2=[];
    disp('Average Reject data') 
    h = waitbar(0,'Average Reject data...');
    for cpt_set=1:1:enum.nset
        for cpt_slc=1:1:enum.datasize(cpt_set).slc
         for cpt_b=1:1:enum.datasize(cpt_set).b     
           for cpt_dir=1:1: enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).nb_dir    
            Mask_Avg=[];   
            if cpt_b>1

                tmpADC=log(Dcm(:,:,cpt_slc,cpt_b,cpt_dir,1:enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg,cpt_set )./repmat(Dcm(:,:,cpt_slc,1,1,1,cpt_set),1,1,1,1,1,enum.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg))/-enum.b(cpt_b); % tmpADC= log(S/S0)/-b-value
                tmpADC(tmpADC>ADC_limit|tmpADC<0.000001)=0;    
                tmpADC(tmpADC>0.000001)=1;

                Mask_Avg=sum(tmpADC,6);

                Dcm2(:,:,cpt_slc,cpt_b,cpt_dir,1,cpt_set)=sum(Dcm(:,:,cpt_slc,cpt_b,cpt_dir,1:enum.dataset(1).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg,cpt_set).*tmpADC,6)./Mask_Avg;
            else
                Dcm2(:,:,cpt_slc,cpt_b,cpt_dir,1,cpt_set)=mean(Dcm(:,:,cpt_slc,cpt_b,cpt_dir,1:enum.dataset(1).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg,cpt_set),6);
            end
            enum2.dataset(cpt_set).slc(cpt_slc).b(cpt_b).dir(cpt_dir).nb_avg=1;
            
           end
         end   
         waitbar(cpt_slc/enum.datasize(cpt_set).slc,h);
        end
    end
   
    Dcm2(isnan(Dcm2))=0;
    Dcm2(isinf(Dcm2))=0;
    close(h);    

end