
clear all
close all
warning off;

gif_mode = false;
rigid_mode = false;
Nrigid_mode = false;
pca_mode = false;
avg_mode = false;
avg2_mode = true;
tmip_mode = false;
inter_mode = false;
trace_mode = true;
mask_mode = true;
ADC_mode = true;
roi_mode = false;
DTI_mode = true;


%%%%%%%%%%%%%%% File Managmeent %%%%%%%%%%%%%%%%%%%%%%
basepath = "/run/media/dogiermann/7290070e-5710-46ab-abc4-09d39f79738b/luigi/DENSE/healthy/";
all_dirs = [
    strcat(basepath,"S10/Ennis_Tr1242_Exvivo/resolve_15seg_B1000_10x10x10/"),
    strcat(basepath,"S10/Ennis_Tr1242_Exvivo/resolve_15seg_B1000_10x10x10_32/"),
    strcat(basepath,"S10/Ennis_Tr1242_Exvivo/resolve_15seg_B1000_10x10x10_33/"),
    strcat(basepath,"S10/Ennis_Tr1242_Exvivo/resolve_15seg_B1000_10x10x10_34/"),
    strcat(basepath,"S10/Ennis_Tr1242_Exvivo/resolve_15seg_B1000_10x10x10_35/"),
    strcat(basepath,"S11/Ennis_Tr1243_Exvivo/resolve_15seg_B1000_10x10x10/"),
    strcat(basepath,"S11/Ennis_Tr1243_Exvivo/resolve_15seg_B1000_10x10x10_31/"),
    strcat(basepath,"S11/Ennis_Tr1243_Exvivo/resolve_15seg_B1000_10x10x10_32/"),
    strcat(basepath,"S11/Ennis_Tr1243_Exvivo/resolve_15seg_B1000_10x10x10_33/"),
    strcat(basepath,"S11/Ennis_Tr1243_Exvivo/resolve_15seg_B1000_10x10x10_34/"),
    strcat(basepath,"S11/Ennis_Tr1243_Exvivo/resolve_15seg_B1000_10x10x10_35/"),
    strcat(basepath,"S12/Ennis_Tr1264_Exvivo/resolve_15seg_B1000_10x10x10/"),
    strcat(basepath,"S12/Ennis_Tr1264_Exvivo/resolve_15seg_B1000_10x10x10_31/"),
    strcat(basepath,"S12/Ennis_Tr1264_Exvivo/resolve_15seg_B1000_10x10x10_32/"),
    strcat(basepath,"S12/Ennis_Tr1264_Exvivo/resolve_15seg_B1000_10x10x10_33/"),
    strcat(basepath,"S12/Ennis_Tr1264_Exvivo/resolve_15seg_B1000_10x10x10_34/"),
    strcat(basepath,"S12/Ennis_Tr1264_Exvivo/resolve_15seg_B1000_10x10x10_35/"),
    strcat(basepath,"S13/Ennis_Tr1241_Exvivo/resolve_15seg_B1000_10x10x10/"),
    strcat(basepath,"S13/Ennis_Tr1241_Exvivo/resolve_15seg_B1000_10x10x10_31/"),
    strcat(basepath,"S13/Ennis_Tr1241_Exvivo/resolve_15seg_B1000_10x10x10_32/"),
    strcat(basepath,"S13/Ennis_Tr1241_Exvivo/resolve_15seg_B1000_10x10x10_33/"),
    strcat(basepath,"S13/Ennis_Tr1241_Exvivo/resolve_15seg_B1000_10x10x10_34/"),
    strcat(basepath,"S13/Ennis_Tr1241_Exvivo/resolve_15seg_B1000_10x10x10_35/"),
    strcat(basepath,"S14/TR1309_Exvivo/resolve1x1x1/"),
    strcat(basepath,"S14/TR1309_Exvivo/resolve1x1x1_15/"),
    strcat(basepath,"S14/TR1309_Exvivo/resolve1x1x1_16/"),
    strcat(basepath,"S14/TR1309_Exvivo/resolve1x1x1_17/"),
    strcat(basepath,"S14/TR1309_Exvivo/resolve1x1x1_18/"),
    strcat(basepath,"S14/TR1309_Exvivo/resolve1x1x1_44/"),
    strcat(basepath,"S15/TR1310_ExVivo/resolve_11seg_B1000_10x10x10/"),
    strcat(basepath,"S15/TR1310_ExVivo/resolve_11seg_B1000_10x10x10_7/"),
    strcat(basepath,"S15/TR1310_ExVivo/resolve_11seg_B1000_10x10x10_8/"),
    strcat(basepath,"S15/TR1310_ExVivo/resolve_11seg_B1000_10x10x10_9/"),
    strcat(basepath,"S15/TR1310_ExVivo/resolve_11seg_B1000_10x10x10_10/"),
    strcat(basepath,"S15/TR1310_ExVivo/resolve_11seg_B1000_10x10x10_11/"),
    strcat(basepath,"S16/TR1311_ExVivo/resolve_11seg_B1000_10x10x10/"),
    strcat(basepath,"S16/TR1311_ExVivo/resolve_11seg_B1000_10x10x10_7/"),
    strcat(basepath,"S16/TR1311_ExVivo/resolve_11seg_B1000_10x10x10_8/"),
    strcat(basepath,"S16/TR1311_ExVivo/resolve_11seg_B1000_10x10x10_9/"),
    strcat(basepath,"S16/TR1311_ExVivo/resolve_11seg_B1000_10x10x10_10/"),
    strcat(basepath,"S16/TR1311_ExVivo/resolve_11seg_B1000_10x10x10_11/"),
    strcat(basepath,"S17/TR1321_ExVivo/Ex_Vivotr1321_Ennis_19000101/resolve_11seg_B1000_10x10x10/"),
    strcat(basepath,"S17/TR1321_ExVivo/Ex_Vivotr1321_Ennis_19000101/resolve_11seg_B1000_10x10x10_8_010902/"),
    strcat(basepath,"S17/TR1321_ExVivo/Ex_Vivotr1321_Ennis_19000101/resolve_11seg_B1000_10x10x10_9_014706/"),
    strcat(basepath,"S17/TR1321_ExVivo/Ex_Vivotr1321_Ennis_19000101/resolve_11seg_B1000_10x10x10_10_035908/"),
    strcat(basepath,"S17/TR1321_ExVivo/Ex_Vivotr1321_Ennis_19000101/resolve_11seg_B1000_10x10x10_11_061139/"),
    strcat(basepath,"S17/TR1321_ExVivo/Ex_Vivotr1321_Ennis_19000101/resolve_11seg_B1000_10x10x10_12_070919/"),
];

for dirid = 1:length(all_dirs)
    dcm_dir = string(all_dirs(dirid));

    cd(dcm_dir);
    mkdir(strcat(dcm_dir,'/Maps'))
    if gif_mode
        mkdir(strcat(dcm_dir,'/Gif'))
    end
    listing = dir(dcm_dir);
    
    %%
    %%%%%%%%%%%%%%% Create Enum and Vol %%%%%%%%%%%%%%%%
    [Dcm enum]= AnalyseDataSet_KM(listing);
    %[Dcm enum]= AnalyseDataSet_forced_KM(listing, [5],[0 350],[1 6],[5 5]);
    enum.dcm_dir=dcm_dir;
    enum.nset=1;
    save(strcat(dcm_dir,'/Maps/RAW.mat'),'Dcm','enum');
    if gif_mode
        Gif_KM(Dcm, enum, 'Raw')
    end
    
    if (enum.dataset.slc(1).b(1).dir(1).nb_avg>enum.dataset.slc(1).b(2).dir(1).nb_avg) % there is more b0 than b-values therefore it's T2 values
        enum.dataset.slc(1).b(1).dir(1).nb_avg=enum.dataset.slc(1).b(2).dir(1).nb_avg;
        DcmB0_T2=Dcm(:,:,1,1,1,enum.dataset.slc(1).b(2).dir(1).nb_avg:end);
        Dcm(:,:,1,1,1,enum.dataset.slc(1).b(2).dir(1).nb_avg:end)=nan;
    end
    
    %%
    %%%%%%%%%%%%%%% Unmosaic %%%%%%%%%
    if enum.mosa>1
        [Dcm, enum]= Demosa_KM(Dcm, enum);
        save(strcat(dcm_dir,'/Maps/Demosa.mat'),'Dcm','enum');
        if gif_mode
            Gif_KM(Dcm, enum, 'Unmosaic')
        end
    end
    
    %%
    %%%%%%%%%%%%%%% Registration Before %%%%%%%%%
    if rigid_mode
        [Dcm]= RigidRegistration_KM2(Dcm, enum);
        save(strcat(dcm_dir,'/Maps/Rigid2.mat'),'Dcm','enum');
        if gif_mode
            Gif_KM(Dcm, enum, 'RigidReg')
        end
    end
    
    if Nrigid_mode
        [Dcm]= NonRigidRegistration_KM(Dcm, enum);
        save(strcat(dcm_dir,'/Maps/NonRigid.mat'),'Dcm','enum');
        if gif_mode
            Gif_KM(Dcm, enum, 'NonRigidReg')
        end
    end
    
    
    
    %%
    %%%%%%%%%%%%%%% PCA %%%%%%%%%
    if pca_mode
        [Dcm ]= VPCA_KM(Dcm,enum,80);
        save(strcat(dcm_dir,'/Maps/PCA.mat'),'Dcm','enum');
    end
    
    %%
    %%%%%%%%%%%%%%% tMIP %%%%%%%%%
    if tmip_mode
        [Dcm enum]= tMIP_KM(Dcm, enum);
        save(strcat(dcm_dir,'/Maps/tMIP.mat'),'Dcm','enum');
         if gif_mode
            Gif_KM(Dcm, enum, 'tMIP')
        end
    end
    
    %%
    %%%%%%%%%%%%%%%% Average %%%%%%%%%
    if avg_mode   
        [Dcm enum]= Average_KM(Dcm, enum); 
        save(strcat(dcm_dir,'/Maps/Average.mat'),'Dcm','enum');
        if gif_mode
            Gif_KM(Dcm, enum, 'Average')
        end
    end
    
    %%
    %%%%%%%%%%%%%%%% Average and reject %%%%%%%%%
    if avg2_mode 
        [Dcm  enum]= Average_and_Reject_KM(Dcm, enum,4e-3);
        save(strcat(dcm_dir,'/Maps/Average_Reject.mat'),'Dcm','enum');
         if gif_mode
            Gif_KM(Dcm, enum, 'Average_Reject')
        end
    end
    %%
    %%%%%%%%%%%%%%% Registration After %%%%%%%%%
    if rigid_mode
        [Dcm]= RigidRegistration_KM(Dcm, enum);
        save(strcat(dcm_dir,'/Maps/RigidAfter.mat'),'Dcm','enum');
        if gif_mode
            Gif_KM(Dcm, enum, 'RigidReg')
        end
    end
    
    %%
    %%%%%%%%%%%%%%% Zero filling interpolation %%%%%%%%%
    if inter_mode
       % [Dcm2 enum]= Depolation_KM(Dcm, enum);
    
       [Dcm enum]= Interpolation_KM(Dcm, enum);
        save(strcat(dcm_dir,'/Maps/Interpolation.mat'),'Dcm','enum');
         if gif_mode
            Gif_KM(Dcm, enum, 'Interpolation')
        end
    end 
    %[Dcm enum]= Mean_KM(Dcm, enum);
    %%
    %%%%%%%%%%%%%%% Create Trace %%%%%%%%%
    if trace_mode
        [Trace enum]= Trace_KM(Dcm, enum,1);
        [Trace_Norm]= Norm_KM(Trace, enum);  
        save(strcat(dcm_dir,'/Maps/Trace.mat'),'Trace','Trace_Norm','enum');  
        if gif_mode
            Gif_KM(Trace, enum, 'Trace')
        end
    end
    
    %%
    %%%%%%%%%%%%%%% Calculate ADC %%%%%%%%%
    if ADC_mode   
        [ADC]= ADC_KM(Trace, enum);
        if enum.datasize.b>2
           ADC(:,:,:,3)=log(Trace(:,:,:,3)./Trace(:,:,:,2)) ./ (enum.b(2)-enum.b(3));
        end
        save(strcat(dcm_dir,'/Maps/ADC.mat'),'ADC');
    end
    
    
    %%
    %%%%%%%%%%%%%%% Create Mask %%%%%%%%%
    if mask_mode
        [Mask]= Mask_KM(Trace(:,:,:,2),60,60000);
        Mask(Mask>0)=1;
        %Dcm=Apply_Mask_KM(Dcm,Mask);
        save(strcat(dcm_dir,'/Maps/Mask.mat'),'Mask','Dcm');  
    end
    %%
    %%%%%%%%%%%%%%% Create ROI %%%%%%%%%
    if roi_mode    
         if isfile(strcat(dcm_dir,'/Maps/ROI.mat'))
         load (strcat(dcm_dir,'/Maps/ROI.mat'));
         %else   
         %    [P_Endo,P_Epi,LV_Mask,Mask_Depth]= ROI_KM(Trace(:,:,:,2));
         %   [Mask_AHA] = ROI2AHA_KM (Dcm, P_Endo, P_Epi);
            %Dcm=Apply_Mask_KM(Dcm,LV_mask);
         %   save(strcat(dcm_dir,'/Maps/ROI.mat'),'P_Endo','P_Epi','LV_Mask','Mask_AHA','Mask_Depth');
         end
    end
    %%
    %%%%%%%%%%%%%%% Calculate Tensor %%%%%%%%%
    if DTI_mode
        [Tensor,EigValue,EigVector,MD,FA,Trace_DTI] = Calc_Tensor_KM(Dcm, enum);
        save(strcat(dcm_dir,'/Maps/DTI.mat'),'Tensor','EigValue','EigVector','MD','FA','Trace_DTI');  
        
       %%%%%%%%%%%%%% Extract HA %%%%%%%%%%%%%
        if roi_mode
                    
            EigVect1=[];
            EigVect2=[];
            EigVect3=[];
            Elevation=[];
            EigVect1(:,:,1:size(EigVector,3),:)=squeeze(EigVector(:,:,:,1,:,1));
            EigVect2(:,:,1:size(EigVector,3),:)=squeeze(EigVector(:,:,:,1,:,2));
            EigVect3(:,:,1:size(EigVector,3),:)=squeeze(EigVector(:,:,:,1,:,3));
            Elevation(:,:,1:size(EigVector,3),:)=squeeze(EigVector(:,:,:,1,3,1));
            
            
            
            %[HA TRA]= HA_KM( EigVect1, Dcm, P_Epi, P_Endo );
            [HA2 TRA2 E2A ]= HA_E2A_KM(EigVect1, EigVect2, LV_Mask, P_Epi, P_Endo);
            %[HA_filter]= HA_Filter_KM(HA,LV_Mask ,Mask_Depth,0);n
    
            [HA_filter2]= HA_Filter_KM(HA2,LV_Mask ,Mask_Depth,0);
            [HA_filter22]= HA_Filter_KM(HA_filter2,LV_Mask ,Mask_Depth,0);
            save(strcat(dcm_dir,'/Maps/HA2.mat'),'EigVect1','EigVect2','EigVect3','Elevation','HA2','TRA2','HA_filter2','E2A');        
            
            % DTI2VTK_KM(EigVector,LV_Mask, enum,1,[],'test',[1 1 1]); % Uncomment to export fiber to
            % VTK
        end
    end
    %%
    % if ADC_mode  && trace_mode 
    %     [Folder]= Recreate_Dicom_Maps_KM(ADC*1e6.*Mask,enum,[],'ADCMap',1015);
    % 	Recreate_Dicom_Maps_KM(Trace*1e1.*Mask,enum,Folder,'TraceMap',1016);
    % end
end

%%
warning on;