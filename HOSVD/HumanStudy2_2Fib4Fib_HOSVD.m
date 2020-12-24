%Analysis of Human Study 2 - Two fibre and Four fibre Blend
%Data: 02-12-2020
%Authors: Chandni Desai, Omar Delannoy-Bruno

%Purpose
%This script will be used to analyze:
    %3. Human Study 2 : Two and Four Fibre blend
    
 %These datasets will be analyzed at the level of 
      %1. ASVs
      %2. mcSEED subsystems
      %3. CAZymes

    
    %% Human Study 2 (workspace: MDZHS4_ASV_HOSVD_analysis.mat)
    
    %Step 1: Concatenate data from days 09,11,25,35,49 into a single three-dimensional object 
    
    %Normalization procedure
        %All data is transformed by log2FC relative to Day 09
        
        %This requires the addition of pseudocounts to the data. 
        %Pseudocount of 0.5 chosen
       
    pseudocount = 0.5;
    HS4data.day09 = Day09 + pseudocount;    %unsupp HiSFLoFV
    HS4data.day11 = Day11 + pseudocount;    %unsupp HiSFLoFV
    HS4data.day25 = Day25 + pseudocount;    %2fibreblend
    HS4data.day35 = Day35 + pseudocount;    %unsupp HiSFLoFV
    HS4data.day49 = Day49 + pseudocount;    %4fibreblend
  
    %Compute log2FC
    
    HS4data.day09_2FC = log2(HS4data.day09./HS4data.day09);
    HS4data.day11_2FC = log2(HS4data.day11./HS4data.day09);
    HS4data.day25_2FC = log2(HS4data.day25./HS4data.day09);
    HS4data.day35_2FC = log2(HS4data.day35./HS4data.day09);
    HS4data.day49_2FC = log2(HS4data.day49./HS4data.day09);

  
    %% ASV : 2 fibre blend
    
    %Concatenate log2FC data for subset into a three-dimensional matrix
    
    HS4data.matrix_2fibreblend(:,:,1) = HS4data.day09_2FC;
    HS4data.matrix_2fibreblend(:,:,2) = HS4data.day11_2FC;
    HS4data.matrix_2fibreblend(:,:,3) = HS4data.day25_2FC;
    HS4data.matrix_2fibreblend(:,:,4) = HS4data.day35_2FC;
  
  
    %HOSVD
    [HS4data.M1_2fibreblend,HS4data.out_2fibreblend,HS4data.thrsh_2fibreblend,HS4data.random_component_2fibreblend] = compute_M1_w_RMT(HS4data.matrix_2fibreblend);
    
    %YOU HAVE NOW COMPLETED HOSVD
    
    %Plot
    for i=1:HS4data.random_component_2fibreblend;
        figure; 
        bar(HS4data.M1_2fibreblend.U{3}(:,i));
    end;
   
    figure; scatter3(HS4data.M1_2fibreblend.U{3}(:,1),HS4data.M1_2fibreblend.U{3}(:,2),HS4data.M1_2fibreblend.U{3}(:,3));
    
    

 %% ASV : 4 fibre blend   
     
    %Concatenate log2FC data for subset into a three-dimensional matrix
    
    HS4data.matrix_4fibreblend(:,:,1) = HS4data.day09_2FC;
    HS4data.matrix_4fibreblend(:,:,2) = HS4data.day11_2FC;
    HS4data.matrix_4fibreblend(:,:,3) = HS4data.day35_2FC;
    HS4data.matrix_4fibreblend(:,:,4) = HS4data.day49_2FC;
  
  
    %HOSVD
    [HS4data.M1_4fibreblend,HS4data.out_4fibreblend,HS4data.thrsh_4fibreblend,HS4data.random_component_4fibreblend] = compute_M1_w_RMT(HS4data.matrix_4fibreblend);
    
    %YOU HAVE NOW COMPLETED HOSVD
    
    %Plot
    for i=1:HS4data.random_component_4fibreblend;
        figure; 
        bar(HS4data.M1_4fibreblend.U{3}(:,i));
    end;
    
    %Result: 
    
    figure; scatter3(HS4data.M1_4fibreblend.U{3}(:,1),HS4data.M1_4fibreblend.U{3}(:,2),HS4data.M1_4fibreblend.U{3}(:,3));
    


    %% 2. mcSEED for Human Study 4 (workspace: MDZHS4_mcSEED_HOSVD.mat)

    %Step 1: Concatenate data from days 09,11,25,35,49 into a single three-dimensional object 
    
    %Normalization procedure
        %All data is transformed by log2FC relative to Day 09
        
        %This requires the addition of pseudocounts to the data. 
        %Pseudocount of 0.5 chosen
        
    pseudocount = 0.5;
    HS4mcseeddata.day09 = Day09 + pseudocount;    %unsupp HiSFLoFV
    HS4mcseeddata.day11 = Day11 + pseudocount;    %unsupp HiSFLoFV
    HS4mcseeddata.day25 = Day25 + pseudocount;    %2fibreblend
    HS4mcseeddata.day35 = Day35 + pseudocount;    %unsupp HiSFLoFV
    HS4mcseeddata.day49 = Day49 + pseudocount;    %4fibreblend
  
    %Compute log2FC
    HS4mcseeddata.day09_2FC = log2(HS4mcseeddata.day09./HS4mcseeddata.day09);
    HS4mcseeddata.day11_2FC = log2(HS4mcseeddata.day11./HS4mcseeddata.day09);
    HS4mcseeddata.day25_2FC = log2(HS4mcseeddata.day25./HS4mcseeddata.day09);
    HS4mcseeddata.day35_2FC = log2(HS4mcseeddata.day35./HS4mcseeddata.day09);
    HS4mcseeddata.day49_2FC = log2(HS4mcseeddata.day49./HS4mcseeddata.day09);
    

    %% mcSEED subsystems : 2 fibre blend 

    %Concatenate log2FC data for subset into a three-dimensional matrix
    HS4mcseeddata.matrix_2fibreblend(:,:,1) = HS4mcseeddata.day09_2FC;
    HS4mcseeddata.matrix_2fibreblend(:,:,2) = HS4mcseeddata.day11_2FC;
    HS4mcseeddata.matrix_2fibreblend(:,:,3) = HS4mcseeddata.day25_2FC;
    HS4mcseeddata.matrix_2fibreblend(:,:,4) = HS4mcseeddata.day35_2FC;
  
  
    %HOSVD on 2 fibre blend
    [HS4mcseeddata.M1_2fibreblend,HS4mcseeddata.out_2fibreblend,HS4mcseeddata.thrsh_2fibreblend,HS4mcseeddata.random_component_2fibreblend] = compute_M1_w_RMT(HS4mcseeddata.matrix_2fibreblend);
    
    %YOU HAVE NOW COMPLETED HOSVD
    
    %Plot
    for i=1:HS4mcseeddata.random_component_2fibreblend;
        figure; 
        bar(HS4mcseeddata.M1_2fibreblend.U{3}(:,i));
    end;
    
    %Result: 
    
    figure; scatter3(HS4mcseeddata.M1_2fibreblend.U{3}(:,1),HS4mcseeddata.M1_2fibreblend.U{3}(:,2),HS4mcseeddata.M1_2fibreblend.U{3}(:,3));
   
    
     %% mcSEED : 4 fibre blend 

    %Concatenate log2FC data for subset into a three-dimensional matrix    
    HS4mcseeddata.matrix_4fibreblend(:,:,1) = HS4mcseeddata.day09_2FC;
    HS4mcseeddata.matrix_4fibreblend(:,:,2) = HS4mcseeddata.day11_2FC;
    HS4mcseeddata.matrix_4fibreblend(:,:,3) = HS4mcseeddata.day35_2FC;
    HS4mcseeddata.matrix_4fibreblend(:,:,4) = HS4mcseeddata.day49_2FC;
  
  
    %HOSVD on 4 fibre blend
    [HS4mcseeddata.M1_4fibreblend,HS4mcseeddata.out_4fibreblend,HS4mcseeddata.thrsh_4fibreblend,HS4mcseeddata.random_component_4fibreblend] = compute_M1_w_RMT(HS4mcseeddata.matrix_4fibreblend);
    
    %YOU HAVE NOW COMPLETED HOSVD
    
    %Plot
    for i=1:HS4mcseeddata.random_component_4fibreblend;
        figure; 
        bar(HS4mcseeddata.M1_4fibreblend.U{3}(:,i));
    end;
    
    %Result: 
    
    figure; scatter3(HS4mcseeddata.M1_4fibreblend.U{3}(:,1),HS4mcseeddata.M1_4fibreblend.U{3}(:,2),HS4mcseeddata.M1_4fibreblend.U{3}(:,3));
   
    


    %% 3. CAZyme Analysis (workspace: MDZHS4_CAZyme_HOSVD_analysis.mat)

    %Step 1: Concatenate data from days 09,11,25,35,49 into a single three-dimensional object 
    
    %Normalization procedure
        %All data is transformed by log2FC relative to Day 09
        
        %This requires the addition of pseudocounts to the data. 
        %Pseudocount of 0.5 chosen
        
        %Rank is the minimum dimension of the 3-D matrix to be constructed
        
    pseudocount = 0.5;
    HS4data.day09 = Day09 + pseudocount;    %unsupp HiSFLoFV
    HS4data.day11 = Day11 + pseudocount;    %unsupp HiSFLoFV
    HS4data.day25 = Day25 + pseudocount;    %2fibreblend
    HS4data.day35 = Day35 + pseudocount;    %unsupp HiSFLoFV
    HS4data.day49 = Day49 + pseudocount;    %4fibreblend
  
    %Compute log2FC
    
    HS4data.day09_2FC = log2(HS4data.day09./HS4data.day09);
    HS4data.day11_2FC = log2(HS4data.day11./HS4data.day09);
    HS4data.day25_2FC = log2(HS4data.day25./HS4data.day09);
    HS4data.day35_2FC = log2(HS4data.day35./HS4data.day09);
    HS4data.day49_2FC = log2(HS4data.day49./HS4data.day09);


    %% CAZymes : 2 fibre blend
    
     
    %Concatenate log2FC data for subset into a three-dimensional matrix
    HS4data.matrix_2fibreblend(:,:,1) = HS4data.day09_2FC;
    HS4data.matrix_2fibreblend(:,:,2) = HS4data.day11_2FC;
    HS4data.matrix_2fibreblend(:,:,3) = HS4data.day25_2FC;
    HS4data.matrix_2fibreblend(:,:,4) = HS4data.day35_2FC;
  
  
    %HOSVD
    [HS4data.M1_2fibreblend,HS4data.out_2fibreblend,HS4data.thrsh_2fibreblend,HS4data.random_component_2fibreblend] = compute_M1_w_RMT(HS4data.matrix_2fibreblend);
    
    %YOU HAVE NOW COMPLETED HOSVD
    
    %Plot
    for i=1:HS4data.random_component_2fibreblend;
        figure; 
        bar(HS4data.M1_2fibreblend.U{3}(:,i));
    end;
    
    %Result: 
    
    figure; scatter3(HS4data.M1_2fibreblend.U{3}(:,1),HS4data.M1_2fibreblend.U{3}(:,2),HS4data.M1_2fibreblend.U{3}(:,3));
    
    
    %% CAZymes : 4 fibre blend
    
    %Concatenate log2FC data for subset into a three-dimensional matrix
    HS4data.matrix_4fibreblend(:,:,1) = HS4data.day09_2FC;
    HS4data.matrix_4fibreblend(:,:,2) = HS4data.day11_2FC;
    HS4data.matrix_4fibreblend(:,:,3) = HS4data.day35_2FC;
    HS4data.matrix_4fibreblend(:,:,4) = HS4data.day49_2FC;
  
  
    %HOSVD
    [HS4data.M1_4fibreblend,HS4data.out_4fibreblend,HS4data.thrsh_4fibreblend,HS4data.random_component_4fibreblend] = compute_M1_w_RMT(HS4data.matrix_4fibreblend);
    
    %YOU HAVE NOW COMPLETED HOSVD
    
    %Plot
    for i=1:HS4data.random_component_4fibreblend;
        figure; 
        bar(HS4data.M1_4fibreblend.U{3}(:,i));
    end;
    
    %Result: 
    
    figure; scatter3(HS4data.M1_4fibreblend.U{3}(:,1),HS4data.M1_4fibreblend.U{3}(:,2),HS4data.M1_4fibreblend.U{3}(:,3));

   