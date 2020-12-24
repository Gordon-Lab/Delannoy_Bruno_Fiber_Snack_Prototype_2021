%Analysis of Human Study 1 - Pea fibre
%Data: 11_8_2019
%Authors: Arjun S. Raman, Chandni Desai

%Purpose
%This script will be used to analyze:
    %Human Study 1 : Pea fibre
    
 %These datasets will be analyzed at the level of 
    %1. ASVs
    %2. mcSEED subsystems
    %3. CAZymes
    
    
    %% Human Study 1 (workspace: MDZHS3_ASV_HOSVD_analysis.mat)
    
    %Step 1: Concatenate data from days 14, 15, 16, 17, 18, 19, 20,
                                       %21, 29, 35, 39, 45, 49
    %into a single three-dimensional object 
    
    %Normalization procedure
        %All data is transformed by log2FC relative to Day 14
        
        %This requires the addition of pseudocounts to the data. 
        %Pseudocount of 0.5 chosen
        
        %Rank is the minimum dimension of the 3-D matrix to be constructed
        
    rank = 13; %Confirm    
    pseudocount = 0.5;
    HS3data.day14 = Day14 + pseudocount;
    HS3data.day15 = Day15 + pseudocount;
    HS3data.day16 = Day16 + pseudocount;
    HS3data.day17 = Day17 + pseudocount;
    HS3data.day18 = Day18 + pseudocount;
    HS3data.day19 = Day19 + pseudocount;
    HS3data.day20 = Day20 + pseudocount;
    HS3data.day21 = Day21 + pseudocount;
    HS3data.day29 = Day29 + pseudocount;
    HS3data.day35 = Day35 + pseudocount;
    HS3data.day39 = Day39 + pseudocount;
    HS3data.day45 = Day45 + pseudocount;
    HS3data.day49 = Day49 + pseudocount;
    
    %Compute log2FC
    
    HS3data.day14_2FC = log2(HS3data.day14./HS3data.day14);
    HS3data.day15_2FC = log2(HS3data.day15./HS3data.day14);
    HS3data.day16_2FC = log2(HS3data.day16./HS3data.day14);
    HS3data.day17_2FC = log2(HS3data.day17./HS3data.day14);
    HS3data.day18_2FC = log2(HS3data.day18./HS3data.day14);
    HS3data.day19_2FC = log2(HS3data.day19./HS3data.day14);
    HS3data.day20_2FC = log2(HS3data.day20./HS3data.day14);
    HS3data.day21_2FC = log2(HS3data.day21./HS3data.day14);
    HS3data.day29_2FC = log2(HS3data.day29./HS3data.day14);
    HS3data.day35_2FC = log2(HS3data.day35./HS3data.day14);
    HS3data.day39_2FC = log2(HS3data.day39./HS3data.day14);
    HS3data.day45_2FC = log2(HS3data.day45./HS3data.day14);
    HS3data.day49_2FC = log2(HS3data.day49./HS3data.day14);
    
    
    %Concatenate all log2FC data into a three-dimensional matrix
    
    HS3data.matrix_full(:,:,1) = HS3data.day14_2FC;
    HS3data.matrix_full(:,:,2) = HS3data.day15_2FC;
    HS3data.matrix_full(:,:,3) = HS3data.day16_2FC;
    HS3data.matrix_full(:,:,4) = HS3data.day17_2FC;
    HS3data.matrix_full(:,:,5) = HS3data.day18_2FC;
    HS3data.matrix_full(:,:,6) = HS3data.day19_2FC;
    HS3data.matrix_full(:,:,7) = HS3data.day20_2FC;
    HS3data.matrix_full(:,:,8) = HS3data.day21_2FC;
    HS3data.matrix_full(:,:,9) = HS3data.day29_2FC;
    HS3data.matrix_full(:,:,10) = HS3data.day35_2FC;
    HS3data.matrix_full(:,:,11) = HS3data.day39_2FC;
    HS3data.matrix_full(:,:,12) = HS3data.day45_2FC;
    HS3data.matrix_full(:,:,13) = HS3data.day49_2FC;
    
    %HOSVD on matrix_full
    [HS3data.M1,HS3data.out,HS3data.thrsh,HS3data.random_component] = compute_M1_w_RMT(HS3data.matrix_full);
    
    %YOU HAVE NOW COMPLETED HOSVD
    
    %Plot
    for i=1:HS3data.random_component;
        figure; 
        bar(HS3data.M1.U{3}(:,i));
    end;
    
    figure; scatter3(HS3data.M1.U{3}(:,1),HS3data.M1.U{3}(:,2),HS3data.M1.U{3}(:,3)); %saved figure
    
    %What are the taxa that determine the variation along each axis?
    
    figure; histogram(HS3data.M1.U{2}(:,2)); ylim([0 10]);
    figure; histogram(HS3data.M1.U{2}(:,3)); ylim([0 10]);
    figure; histogram(HS3data.M1.U{2}(:,4)); ylim([0 10]);
    
    HS3data.taxa_projections.TC2 = HS3data.M1.U{2}(:,2);
    HS3data.taxa_projections.TC3 = HS3data.M1.U{2}(:,3);
    HS3data.taxa_projections.TC4 = HS3data.M1.U{2}(:,4);
    

    
    %% 2. mcSEED Analysis (workspace: MDZHS3_mcSEED_HOSVD.mat)
    
    %Step 1: Concatenate data from days 14, 15, 16, 17, 18, 19, 20,
                                       %21, 29, 35, 39, 45, 49
    %into a single three-dimensional object 
    
    %Normalization procedure
        %All data is transformed by log2FC relative to Day 14
        
        %This requires the addition of pseudocounts to the data. 
        %Pseudocount of 0.5 chosen
        
        %Rank is the minimum dimension of the 3-D matrix to be constructed
        
    rank = 13; %Confirm    
    pseudocount = 0.5;
    HS3mcseed.day14 = Day14 + pseudocount;
    HS3mcseed.day15 = Day15 + pseudocount;
    HS3mcseed.day16 = Day16 + pseudocount;
    HS3mcseed.day17 = Day17 + pseudocount;
    HS3mcseed.day18 = Day18 + pseudocount;
    HS3mcseed.day19 = Day19 + pseudocount;
    HS3mcseed.day20 = Day20 + pseudocount;
    HS3mcseed.day21 = Day21 + pseudocount;
    HS3mcseed.day29 = Day29 + pseudocount;
    HS3mcseed.day35 = Day35 + pseudocount;
    HS3mcseed.day39 = Day39 + pseudocount;
    HS3mcseed.day45 = Day45 + pseudocount;
    HS3mcseed.day49 = Day49 + pseudocount;
    
    %Compute log2FC
    
    HS3mcseed.day14_2FC = log2(HS3mcseed.day14./HS3mcseed.day14);
    HS3mcseed.day15_2FC = log2(HS3mcseed.day15./HS3mcseed.day14);
    HS3mcseed.day16_2FC = log2(HS3mcseed.day16./HS3mcseed.day14);
    HS3mcseed.day17_2FC = log2(HS3mcseed.day17./HS3mcseed.day14);
    HS3mcseed.day18_2FC = log2(HS3mcseed.day18./HS3mcseed.day14);
    HS3mcseed.day19_2FC = log2(HS3mcseed.day19./HS3mcseed.day14);
    HS3mcseed.day20_2FC = log2(HS3mcseed.day20./HS3mcseed.day14);
    HS3mcseed.day21_2FC = log2(HS3mcseed.day21./HS3mcseed.day14);
    HS3mcseed.day29_2FC = log2(HS3mcseed.day29./HS3mcseed.day14);
    HS3mcseed.day35_2FC = log2(HS3mcseed.day35./HS3mcseed.day14);
    HS3mcseed.day39_2FC = log2(HS3mcseed.day39./HS3mcseed.day14);
    HS3mcseed.day45_2FC = log2(HS3mcseed.day45./HS3mcseed.day14);
    HS3mcseed.day49_2FC = log2(HS3mcseed.day49./HS3mcseed.day14);
    
    
    %Concatenate all log2FC data into a three-dimensional matrix
    
    HS3mcseed.matrix_full(:,:,1) = HS3mcseed.day14_2FC;
    HS3mcseed.matrix_full(:,:,2) = HS3mcseed.day15_2FC;
    HS3mcseed.matrix_full(:,:,3) = HS3mcseed.day16_2FC;
    HS3mcseed.matrix_full(:,:,4) = HS3mcseed.day17_2FC;
    HS3mcseed.matrix_full(:,:,5) = HS3mcseed.day18_2FC;
    HS3mcseed.matrix_full(:,:,6) = HS3mcseed.day19_2FC;
    HS3mcseed.matrix_full(:,:,7) = HS3mcseed.day20_2FC;
    HS3mcseed.matrix_full(:,:,8) = HS3mcseed.day21_2FC;
    HS3mcseed.matrix_full(:,:,9) = HS3mcseed.day29_2FC;
    HS3mcseed.matrix_full(:,:,10) = HS3mcseed.day35_2FC;
    HS3mcseed.matrix_full(:,:,11) = HS3mcseed.day39_2FC;
    HS3mcseed.matrix_full(:,:,12) = HS3mcseed.day45_2FC;
    HS3mcseed.matrix_full(:,:,13) = HS3mcseed.day49_2FC;
    
    %HOSVD on matrix_full
    [HS3mcseed.M1,HS3mcseed.out,HS3mcseed.thrsh,HS3mcseed.random_component] = compute_M1_w_RMT(HS3mcseed.matrix_full);
    
    %YOU HAVE NOW COMPLETED HOSVD

    
    
     %% 3. CAZymes analysis (workspace: MDZHS3_CAZymes_HOSVD.mat)
    
    %Step 1: Concatenate data from days 14, 15, 16, 17, 18, 19, 20,
                                       %21, 29, 35, 39, 45, 49
    %into a single three-dimensional object 
    
    %Normalization procedure
        %All data is transformed by log2FC relative to Day 14
        
        %This requires the addition of pseudocounts to the data. 
        %Pseudocount of 0.5 chosen
        
        %Rank is the minimum dimension of the 3-D matrix to be constructed
        
    rank = 13; %Confirm    
    pseudocount = 0.5;
    HS3cazyme.day14 = Day14 + pseudocount;
    HS3cazyme.day15 = Day15 + pseudocount;
    HS3cazyme.day16 = Day16 + pseudocount;
    HS3cazyme.day17 = Day17 + pseudocount;
    HS3cazyme.day18 = Day18 + pseudocount;
    HS3cazyme.day19 = Day19 + pseudocount;
    HS3cazyme.day20 = Day20 + pseudocount;
    HS3cazyme.day21 = Day21 + pseudocount;
    HS3cazyme.day29 = Day29 + pseudocount;
    HS3cazyme.day35 = Day35 + pseudocount;
    HS3cazyme.day39 = Day39 + pseudocount;
    HS3cazyme.day45 = Day45 + pseudocount;
    HS3cazyme.day49 = Day49 + pseudocount;
    
    %Compute log2FC
    
    HS3cazyme.day14_2FC = log2(HS3cazyme.day14./HS3cazyme.day14);
    HS3cazyme.day15_2FC = log2(HS3cazyme.day15./HS3cazyme.day14);
    HS3cazyme.day16_2FC = log2(HS3cazyme.day16./HS3cazyme.day14);
    HS3cazyme.day17_2FC = log2(HS3cazyme.day17./HS3cazyme.day14);
    HS3cazyme.day18_2FC = log2(HS3cazyme.day18./HS3cazyme.day14);
    HS3cazyme.day19_2FC = log2(HS3cazyme.day19./HS3cazyme.day14);
    HS3cazyme.day20_2FC = log2(HS3cazyme.day20./HS3cazyme.day14);
    HS3cazyme.day21_2FC = log2(HS3cazyme.day21./HS3cazyme.day14);
    HS3cazyme.day29_2FC = log2(HS3cazyme.day29./HS3cazyme.day14);
    HS3cazyme.day35_2FC = log2(HS3cazyme.day35./HS3cazyme.day14);
    HS3cazyme.day39_2FC = log2(HS3cazyme.day39./HS3cazyme.day14);
    HS3cazyme.day45_2FC = log2(HS3cazyme.day45./HS3cazyme.day14);
    HS3cazyme.day49_2FC = log2(HS3cazyme.day49./HS3cazyme.day14);
    
    
    %Concatenate all log2FC data into a three-dimensional matrix
    
    HS3cazyme.matrix_full(:,:,1) = HS3cazyme.day14_2FC;
    HS3cazyme.matrix_full(:,:,2) = HS3cazyme.day15_2FC;
    HS3cazyme.matrix_full(:,:,3) = HS3cazyme.day16_2FC;
    HS3cazyme.matrix_full(:,:,4) = HS3cazyme.day17_2FC;
    HS3cazyme.matrix_full(:,:,5) = HS3cazyme.day18_2FC;
    HS3cazyme.matrix_full(:,:,6) = HS3cazyme.day19_2FC;
    HS3cazyme.matrix_full(:,:,7) = HS3cazyme.day20_2FC;
    HS3cazyme.matrix_full(:,:,8) = HS3cazyme.day21_2FC;
    HS3cazyme.matrix_full(:,:,9) = HS3cazyme.day29_2FC;
    HS3cazyme.matrix_full(:,:,10) = HS3cazyme.day35_2FC;
    HS3cazyme.matrix_full(:,:,11) = HS3cazyme.day39_2FC;
    HS3cazyme.matrix_full(:,:,12) = HS3cazyme.day45_2FC;
    HS3cazyme.matrix_full(:,:,13) = HS3cazyme.day49_2FC;
    
    %HOSVD on matrix_full
    [HS3cazyme.M1,HS3cazyme.out,HS3cazyme.thrsh,HS3cazyme.random_component] = compute_M1_w_RMT(HS3cazyme.matrix_full);
    
    %YOU HAVE NOW COMPLETED HOSVD

    
    %What are the taxa that determine the variation along each axis?
    HS3cazyme.taxa_projections.TC1 = HS3cazyme.M1.U{2}(:,1);
    HS3cazyme.taxa_projections.TC2 = HS3cazyme.M1.U{2}(:,2);
    HS3cazyme.taxa_projections.TC3 = HS3cazyme.M1.U{2}(:,3);
    HS3cazyme.taxa_projections.TC3 = HS3cazyme.M1.U{2}(:,4);
    
    
    
