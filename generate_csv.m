clear all
clc
close all
% % % writeCSVFileNew(output,outputName,neighborMatrix);
matRoot='.//MAT1//';
IMCRoot='.//IMC//';

matPath=dir(matRoot);

stain_file_name={'147Sm_CollagenI.ome.tiff','161Dy_CD20.ome.tiff','174Yb_CD57.ome.tiff',...                                                          
'162Dy_CD11c.ome.tiff','175Lu_C1QC.ome.tiff', '148Nd_b2m.ome.tiff' ,'115In_CCR7.ome.tiff','149Sm_CD31.ome.tiff','163Dy_CD15.ome.tiff', '176Yb_Pan-cytokeratin.ome.tiff',...  
'150Nd_E-cadherin.ome.tiff', '164Dy_Granzyme_B.ome.tiff', '151Eu_CD127.ome.tiff', '165Ho_PD1.ome.tiff',...            
'152Sm_CD56.ome.tiff','166Er_Ki67.ome.tiff','194Pt_aSMA.ome.tiff',...             
'153Eu_CD7.ome.tiff', '167Er_CD1c.ome.tiff', '198Pt_Vimentin.ome.tiff',...         
'141Pr_CD14.ome.tiff', '154Sm_CD163.ome.tiff', '168Er_HLA-DR.ome.tiff',...            
'142Nd_FOXP3.ome.tiff', '155Gd_IDO1.ome.tiff', '169Tm_CD45RA.ome.tiff',...        
'143Nd_CD16.ome.tiff','156Gd_PDL1.ome.tiff','170Er_CD3.ome.tiff','89Y_CD45.ome.tiff',...               
'144Nd_HLAI.ome.tiff', '158Gd_VEGFA.ome.tiff','171Yb_TNFa.ome.tiff',...                
'145Nd_CD4.ome.tiff','159Tb_CD68.ome.tiff', '172Yb_CD27.ome.tiff',...             
'146Nd_CD8.ome.tiff', '160Gd_CD11b.ome.tiff', '173Yb_CD45RO.ome.tiff'};    

for i=3:length(matPath)
    temp=matPath(i).name;
    fileName=strcat(matRoot,temp);
    load(fileName)
    stainValue=zeros(max(max(segmentaion)),length(stain_file_name));
    nuclei_coodinate=zeros(max(max(segmentaion)),2);
    temp_split1=regexp(temp,'_','split');
    sampleName=temp_split1{1};
    if length(temp_split1)==4
       temp1=temp_split1{4};
       roiName=strcat(temp_split1{2},'_',temp_split1{3},'_',temp1(1:length(temp1)-4));
    else
       temp1=temp_split1{5};
       roiName=strcat(temp_split1{2},'_',temp_split1{3},'_',temp_split1{4},'_',temp1(1:length(temp1)-4));
    end   
    IMCPath=strcat(IMCRoot,sampleName,'//',roiName);
    outputName=strcat(sampleName,'_',roiName);
    strcat('.//csvFile//',outputName,'.csv')
    imgPath=dir(IMCPath);
    stain_id=0;
    stain_name={};
    for j=3:length(imgPath)
        
        tiffFileName=imgPath(j).name;
        if  sum(contains(stain_file_name,tiffFileName))==1
% % %             tiffFileName
            stain_id=stain_id+1;
            stain_name{stain_id}=tiffFileName;
            stainIndex=find(contains(stain_file_name,tiffFileName)==1);
            inputImageFileRoot=strcat(IMCPath,'\\',tiffFileName);
            im = tiffread2(inputImageFileRoot);
            dataOrigin=uint8(im.data);
            tk=0;
            for k=1:max(max(segmentaion))
                if length(find(segmentaion==k))>0
                tk=tk+1;
                stainValue(tk,stain_id)=mean(dataOrigin(find(segmentaion==k)));
                [A,B]=find(segmentaion==k);
                nuclei_coodinate(tk,:)=[mean(B),mean(A)];
                end
            end

        end
    end

    output=[nuclei_coodinate,stainValue];
    neighborMatrix=getNucleiCoordinate(nuclei_coodinate);
    
    writeCSVFileNew(output,outputName,neighborMatrix);
end
