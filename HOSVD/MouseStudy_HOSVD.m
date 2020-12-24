%Analysis of Mouse Study
%Date: 11_6_2019
%Authors: Arjun S. Raman, Chandni Desai

%Purpose
%This script will be used to perform HOSVD on:
    %Mouse study (57 mice subject to gavage with 9 donor microbiota
    %subject to engineered dietary conditions)
    
 %Analyzed at the level of 
    %1. ASVs
    %2. mcSEED subsystems
    %3. CAZymes
    
    
    %% Mouse Study (workspace: analysis_mousestudy.mat)
    
    %Step 1: Concatenate data from days 14, 24, 34, 44, 54, 64 into a
    %single three-dimensional object 
    
    %Normalization procedure
        %All data is transformed by log2FC relative to Day 14 (the first
        %day of sampling)
        
        %This requires the addition of pseudocounts to the data. 
        %Pseudocount of 0.5 chosen
        
        %Rank is the minimum dimension of the 3-D matrix to be constructed
        
    rank = 6;    
    pseudocount = 0.5;
    data.day14 = Day14 + pseudocount;
    data.day24 = Day24 + pseudocount;
    data.day34 = Day34 + pseudocount;
    data.day44 = Day44 + pseudocount;
    data.day54 = Day54 + pseudocount;
    data.day64 = Day64 + pseudocount;
    
    %Compute log2FC
    
    data.day14_2FC = log2(data.day14./data.day14);
    data.day24_2FC = log2(data.day24./data.day14);
    data.day34_2FC = log2(data.day34./data.day14);
    data.day44_2FC = log2(data.day44./data.day14);
    data.day54_2FC = log2(data.day54./data.day14);
    data.day64_2FC = log2(data.day64./data.day14);
    
    
    %Concatenate all log2FC data into a three-dimensional matrix
    
    data.matrix_full(:,:,1) = data.day14_2FC;
    data.matrix_full(:,:,2) = data.day24_2FC;
    data.matrix_full(:,:,3) = data.day34_2FC;
    data.matrix_full(:,:,4) = data.day44_2FC;
    data.matrix_full(:,:,5) = data.day54_2FC;
    data.matrix_full(:,:,6) = data.day64_2FC;
    
        
    %% Pea fibre analysis
    
    %Isolate data for pea fibre 
    data.matrix_pf = data.matrix_full(:,:,1:3);
    
    %HOSVD on pea fibre
    [data.M1_peafibre,data.out_peafibre,data.thrsh_peafibre,data.random_component_peafibre] = compute_M1_w_RMT(data.matrix_pf);
    
    %YOU HAVE NOW COMPLETED HOSVD
    
    %Plot TC1 and TC2 through time
    for i=1:data.random_component_peafibre;
        figure; 
        bar(data.M1_peafibre.U{3}(:,i));
    end;
    
    figure; scatter(data.M1_peafibre.U{3}(:,1),data.M1_peafibre.U{3}(:,2));
    
    %What are the taxa that determine the variation along each axis?
    figure; histogram(data.M1_peafibre.U{2}(:,1)); ylim([0 10]);
    figure; histogram(data.M1_peafibre.U{2}(:,2)); ylim([0 10]);

    data.taxa_projections_peafibre.TC1 = data.M1_peafibre.U{2}(:,1);
    data.taxa_projections_peafibre.TC2 = data.M1_peafibre.U{2}(:,2);
   
    
    %% Orange fibre analysis
   
    %Compute log2FC with Day 34 as baseline 
    data.day34_2FC_byd34 = log2(data.day34./data.day34);
    data.day44_2FC_byd34 = log2(data.day44./data.day34);
    data.day54_2FC_byd34 = log2(data.day54./data.day34);
     
    %Isolate data for orange fibre 
    data.matrix_of(:,:,1) = data.day34_2FC_byd34;
    data.matrix_of(:,:,2) = data.day44_2FC_byd34;
    data.matrix_of(:,:,3) = data.day54_2FC_byd34;
   
    %HOSVD on orange fibre
[data.M1_orangefibre,data.out_orangefibre,data.thrsh_orangefibre,data.random_component_orangefibre] = compute_M1_w_RMT(data.matrix_of);
    
    %YOU HAVE NOW COMPLETED HOSVD
    
    %Plot TC1 and TC2 through time
    for i=1:data.random_component_orangefibre;
        figure; 
        bar(data.M1_orangefibre.U{3}(:,i));
    end;
    
    figure; scatter(data.M1_orangefibre.U{3}(:,1),data.M1_orangefibre.U{3}(:,2));
    
    %What are the taxa that determine the variation along each axis?
    figure; histogram(data.M1_orangefibre.U{2}(:,1)); ylim([0 10]);
    figure; histogram(data.M1_orangefibre.U{2}(:,2)); ylim([0 10]);
    
    data.taxa_projections_orangefibre.TC1 = data.M1_orangefibre.U{2}(:,1);
    data.taxa_projections_orangefibre.TC2 = data.M1_orangefibre.U{2}(:,2);
   
      
    %% Barley bran analysis
   
    %Compute log2FC with Day 34 as baseline 
    data.day34_2FC_byd34 = log2(data.day34./data.day34);
    data.day54_2FC_byd34 = log2(data.day54./data.day34);
    data.day64_2FC_byd34 = log2(data.day64./data.day34);
 
    %Isolate data for barley bran 
    data.matrix_bb(:,:,1) = data.day34_2FC_byd34;
    data.matrix_bb(:,:,2) = data.day54_2FC_byd34;
    data.matrix_bb(:,:,3) = data.day64_2FC_byd34;

   
    %HOSVD on barley bran
    [data.M1_barleybran,data.out_barleybran,data.thrsh_barleybran,data.random_component_barleybran] = compute_M1_w_RMT(data.matrix_bb);
    
    %YOU HAVE NOW COMPLETED HOSVD
    
    %Plot TC1 and TC2 through time
    for i=1:data.random_component_barleybran;
        figure; 
        bar(data.M1_barleybran.U{3}(:,i));
    end;
    
    figure; scatter(data.M1_barleybran.U{3}(:,1),data.M1_barleybran.U{3}(:,2));
    
    %What are the taxa that determine the variation along each axis?
    figure; histogram(data.M1_barleybran.U{2}(:,1)); ylim([0 10]);
    figure; histogram(data.M1_barleybran.U{2}(:,2)); ylim([0 10]);
    
    data.taxa_projections_barleybran.TC1 = data.M1_barleybran.U{2}(:,1);
    data.taxa_projections_barleybran.TC2 = data.M1_barleybran.U{2}(:,2);
   
    
  
    %% mcSEED analysis, mice (workspace: analysis_mousestudy_mcSEED.mat)
    
    %Step 1: Concatenate data from days 14, 24, 34, 44, 54, 64 into a
    %single three-dimensional object 
    
    %Normalization procedure
        %All data is transformed by log2FC relative to Day 14 (the first
        %day of sampling)
        
        %This requires the addition of pseudocounts to the data. 
        %Pseudocount of 0.5 chosen
        
        %Rank is the minimum dimension of the 3-D matrix to be constructed
        
    rank = 6;    
    pseudocount = 0.5;
    mcseeddata.day14 = Day14 + pseudocount;
    mcseeddata.day24 = Day24 + pseudocount;
    mcseeddata.day34 = Day34 + pseudocount;
    mcseeddata.day44 = Day44 + pseudocount;
    mcseeddata.day54 = Day54 + pseudocount;
    mcseeddata.day64 = Day64 + pseudocount;
    
    %Compute log2FC
    
    mcseeddata.day14_2FC = log2(mcseeddata.day14./mcseeddata.day14);
    mcseeddata.day24_2FC = log2(mcseeddata.day24./mcseeddata.day14);
    mcseeddata.day34_2FC = log2(mcseeddata.day34./mcseeddata.day14);
    mcseeddata.day44_2FC = log2(mcseeddata.day44./mcseeddata.day14);
    mcseeddata.day54_2FC = log2(mcseeddata.day54./mcseeddata.day14);
    mcseeddata.day64_2FC = log2(mcseeddata.day64./mcseeddata.day14);
    
    
    %Concatenate all log2FC data into a three-dimensional matrix
    
    mcseeddata.matrix_full(:,:,1) = mcseeddata.day14_2FC;
    mcseeddata.matrix_full(:,:,2) = mcseeddata.day24_2FC;
    mcseeddata.matrix_full(:,:,3) = mcseeddata.day34_2FC;
    mcseeddata.matrix_full(:,:,4) = mcseeddata.day44_2FC;
    mcseeddata.matrix_full(:,:,5) = mcseeddata.day54_2FC;
    mcseeddata.matrix_full(:,:,6) = mcseeddata.day64_2FC;
    
    
    %% mcSEED analysis, pea fibre, mice
    mcseeddata.matrix_peafibre = mcseeddata.matrix_full(:,:,1:3);
    
     %HOSVD on pea fibre
    [mcseeddata.M1_peafibre,mcseeddata.out_peafibre,mcseeddata.thrsh_peafibre,mcseeddata.random_component_peafibre] = compute_M1_w_RMT(mcseeddata.matrix_peafibre);
    
    %YOU HAVE NOW COMPLETED HOSVD
    
    %Plot TC1 and TC2 through time
    for i=1:mcseeddata.random_component_peafibre;
        figure; 
        bar(mcseeddata.M1_peafibre.U{3}(:,i));
    end;
    
    figure; scatter(mcseeddata.M1_peafibre.U{3}(:,1),mcseeddata.M1_peafibre.U{3}(:,2));
    
    %What are the taxa that determine the variation along each axis?
    figure; histogram(mcseeddata.M1_peafibre.U{2}(:,1)); ylim([0 10]);
    figure; histogram(mcseeddata.M1_peafibre.U{2}(:,2)); ylim([0 10]);

    mcseeddata.taxa_projections_peafibre.TC1 = mcseeddata.M1_peafibre.U{2}(:,1);
    mcseeddata.taxa_projections_peafibre.TC2 = mcseeddata.M1_peafibre.U{2}(:,2);
   
     
   %% mcSEED analysis, orange fibre, mice
    
    %Compute log2FC with Day 34 as baseline 
    mcseeddata.day34_2FC_byd34 = log2(mcseeddata.day34./mcseeddata.day34);
    mcseeddata.day44_2FC_byd34 = log2(mcseeddata.day44./mcseeddata.day34);
    mcseeddata.day54_2FC_byd34 = log2(mcseeddata.day54./mcseeddata.day34);
     
    %Isolate data for orange fibre 
    mcseeddata.matrix_of(:,:,1) = mcseeddata.day34_2FC_byd34;
    mcseeddata.matrix_of(:,:,2) = mcseeddata.day44_2FC_byd34;
    mcseeddata.matrix_of(:,:,3) = mcseeddata.day54_2FC_byd34;
    
    %HOSVD on orange fibre
[mcseeddata.M1_orangefibre,mcseeddata.out_orangefibre,mcseeddata.thrsh_orangefibre,mcseeddata.random_component_orangefibre] = compute_M1_w_RMT(mcseeddata.matrix_of);
    
    %YOU HAVE NOW COMPLETED HOSVD
    
    %Plot TC1 and TC2 through time
    for i=1:mcseeddata.random_component_orangefibre;
        figure; 
        bar(mcseeddata.M1_orangefibre.U{3}(:,i));
    end;
    
    figure; scatter(mcseeddata.M1_orangefibre.U{3}(:,1),mcseeddata.M1_orangefibre.U{3}(:,2));
    
    %What are the taxa that determine the variation along each axis?
    figure; histogram(mcseeddata.M1_orangefibre.U{2}(:,1)); ylim([0 10]);
    figure; histogram(mcseeddata.M1_orangefibre.U{2}(:,2)); ylim([0 10]);
    
    mcseeddata.taxa_projections_orangefibre.TC1 = mcseeddata.M1_orangefibre.U{2}(:,1);
    mcseeddata.taxa_projections_orangefibre.TC2 = mcseeddata.M1_orangefibre.U{2}(:,2);
   

    %% Barley bran analysis
   
    %Compute log2FC with Day 34 as baseline 
    mcseeddata.day34_2FC_byd34 = log2(mcseeddata.day34./mcseeddata.day34);
    mcseeddata.day54_2FC_byd34 = log2(mcseeddata.day54./mcseeddata.day34);
    mcseeddata.day64_2FC_byd34 = log2(mcseeddata.day64./mcseeddata.day34);
 
    %Isolate mcseeddata for barley bran 
    mcseeddata.matrix_bb(:,:,1) = mcseeddata.day34_2FC_byd34;
    mcseeddata.matrix_bb(:,:,2) = mcseeddata.day54_2FC_byd34;
    mcseeddata.matrix_bb(:,:,3) = mcseeddata.day64_2FC_byd34;
 
   
    %HOSVD on barley bran  [mcseeddata.M1_barleybran,mcseeddata.out_barleybran,mcseeddata.thrsh_barleybran,mcseeddata.random_component_barleybran] = compute_M1_w_RMT(mcseeddata.matrix_bb);
    
    %YOU HAVE NOW COMPLETED HOSVD
    
    %Plot TC1 and TC2 through time
    for i=1:mcseeddata.random_component_barleybran;
        figure; 
        bar(mcseeddata.M1_barleybran.U{3}(:,i));
    end;
  
    figure; scatter(mcseeddata.M1_barleybran.U{3}(:,1),mcseeddata.M1_barleybran.U{3}(:,2));
    
    %What are the taxa that determine the variation along each axis?
    figure; histogram(mcseeddata.M1_barleybran.U{2}(:,1)); ylim([0 10]);
    figure; histogram(mcseeddata.M1_barleybran.U{2}(:,2)); ylim([0 10]);
 
    mcseeddata.taxa_projections_barleybran.TC1 = mcseeddata.M1_barleybran.U{2}(:,1);
    mcseeddata.taxa_projections_barleybran.TC2 = mcseeddata.M1_barleybran.U{2}(:,2);
   

    %% CAZyme analysis, mice
    
    %Step 1: Concatenate data from days 14, 24, 34, 44, 54, 64 into a
    %single three-dimensional object 
    
    %Normalization procedure
        %All data is transformed by log2FC relative to Day 14 (the first
        %day of sampling)
        
        %This requires the addition of pseudocounts to the data. 
        %Pseudocount of 0.5 chosen
        
        %Rank is the minimum dimension of the 3-D matrix to be constructed
        
    rank = 6;    
    pseudocount = 0.5;
    cazymedata.day14 = Day14 + pseudocount;
    cazymedata.day24 = Day24 + pseudocount;
    cazymedata.day34 = Day34 + pseudocount;
    cazymedata.day44 = Day44 + pseudocount;
    cazymedata.day54 = Day54 + pseudocount;
    cazymedata.day64 = Day64 + pseudocount;
    
    %Compute log2FC
    
    cazymedata.day14_2FC = log2(cazymedata.day14./cazymedata.day14);
    cazymedata.day24_2FC = log2(cazymedata.day24./cazymedata.day14);
    cazymedata.day34_2FC = log2(cazymedata.day34./cazymedata.day14);
    cazymedata.day44_2FC = log2(cazymedata.day44./cazymedata.day14);
    cazymedata.day54_2FC = log2(cazymedata.day54./cazymedata.day14);
    cazymedata.day64_2FC = log2(cazymedata.day64./cazymedata.day14);
    
    
    %Concatenate all log2FC data into a three-dimensional matrix
    
    cazymedata.matrix_full(:,:,1) = cazymedata.day14_2FC;
    cazymedata.matrix_full(:,:,2) = cazymedata.day24_2FC;
    cazymedata.matrix_full(:,:,3) = cazymedata.day34_2FC;
    cazymedata.matrix_full(:,:,4) = cazymedata.day44_2FC;
    cazymedata.matrix_full(:,:,5) = cazymedata.day54_2FC;
    cazymedata.matrix_full(:,:,6) = cazymedata.day64_2FC;
    
    
    %% cazyme analysis, pea fibre, mice
    
    cazymedata.matrix_peafibre = cazymedata.matrix_full(:,:,1:3);
    
     %HOSVD on pea fibre
[cazymedata.M1_peafibre,cazymedata.out_peafibre,cazymedata.thrsh_peafibre,cazymedata.random_component_peafibre] = compute_M1_w_RMT(cazymedata.matrix_peafibre);
    
    %YOU HAVE NOW COMPLETED HOSVD
    
    %Plot TC1 and TC2 through time
    for i=1:cazymedata.random_component_peafibre;
        figure; 
        bar(cazymedata.M1_peafibre.U{3}(:,i));
    end;
    
    %Result: Plot a two-dimensional plot
    
    figure; scatter(cazymedata.M1_peafibre.U{3}(:,1),cazymedata.M1_peafibre.U{3}(:,2)); 
    
    %What are the taxa that determine the variation along each axis?
    figure; histogram(cazymedata.M1_peafibre.U{2}(:,1)); ylim([0 10]);
    figure; histogram(cazymedata.M1_peafibre.U{2}(:,2)); ylim([0 10]);
 
    
    cazymedata.taxa_projections_peafibre.TC1 = cazymedata.M1_peafibre.U{2}(:,1);
    cazymedata.taxa_projections_peafibre.TC2 = cazymedata.M1_peafibre.U{2}(:,2);
   

     %% cazyme analysis, orange fibre, mice
      
    %Compute log2FC with Day 34 as baseline 
    cazymedata.day34_2FC_byd34 = log2(cazymedata.day34./cazymedata.day34);
    cazymedata.day44_2FC_byd34 = log2(cazymedata.day44./cazymedata.day34);
    cazymedata.day54_2FC_byd34 = log2(cazymedata.day54./cazymedata.day34);
     
    %Isolate data for orange fibre 
    cazymedata.matrix_of(:,:,1) = cazymedata.day34_2FC_byd34;
    cazymedata.matrix_of(:,:,2) = cazymedata.day44_2FC_byd34;
    cazymedata.matrix_of(:,:,3) = cazymedata.day54_2FC_byd34;
    
     %HOSVD on orange fibre
[cazymedata.M1_orangefibre,cazymedata.out_orangefibre,cazymedata.thrsh_orangefibre,cazymedata.random_component_orangefibre] = compute_M1_w_RMT(cazymedata.matrix_of);
    
    %YOU HAVE NOW COMPLETED HOSVD
    
    %Plot TC1 and TC2 through time
    for i=1:cazymedata.random_component_orangefibre;
        figure; 
        bar(cazymedata.M1_orangefibre.U{3}(:,i));
    end;
    
    %Result: Plot a two-dimensional plot
    
    figure; scatter(cazymedata.M1_orangefibre.U{3}(:,1),cazymedata.M1_orangefibre.U{3}(:,2));
    
    %What are the taxa that determine the variation along each axis?
    figure; histogram(cazymedata.M1_orangefibre.U{2}(:,1)); ylim([0 10]);
    figure; histogram(cazymedata.M1_orangefibre.U{2}(:,2)); ylim([0 10]);
 
    cazymedata.taxa_projections_orangefibre.TC1 = cazymedata.M1_orangefibre.U{2}(:,1);
    cazymedata.taxa_projections_orangefibre.TC2 = cazymedata.M1_orangefibre.U{2}(:,2);
  

    %% Barley bran analysis
   
    %Compute log2FC with Day 34 as baseline 
    cazymedata.day34_2FC_byd34 = log2(cazymedata.day34./cazymedata.day34);
    cazymedata.day54_2FC_byd34 = log2(cazymedata.day54./cazymedata.day34);
    cazymedata.day64_2FC_byd34 = log2(cazymedata.day64./cazymedata.day34);
 
    %Isolate cazymedata for barley bran 
    cazymedata.matrix_bb(:,:,1) = cazymedata.day34_2FC_byd34;
    cazymedata.matrix_bb(:,:,2) = cazymedata.day54_2FC_byd34;
    cazymedata.matrix_bb(:,:,3) = cazymedata.day64_2FC_byd34;
 
   
    %HOSVD on barley bran
[cazymedata.M1_barleybran,cazymedata.out_barleybran,cazymedata.thrsh_barleybran,cazymedata.random_component_barleybran] = compute_M1_w_RMT(cazymedata.matrix_bb);
    
    %YOU HAVE NOW COMPLETED HOSVD
    
    %Plot TC1 and TC2 through time
    for i=1:cazymedata.random_component_barleybran;
        figure; 
        bar(cazymedata.M1_barleybran.U{3}(:,i));
    end;
  
    figure; scatter(cazymedata.M1_barleybran.U{3}(:,1),cazymedata.M1_barleybran.U{3}(:,2));
    
    %What are the taxa that determine the variation along each axis?
    
    figure; histogram(cazymedata.M1_barleybran.U{2}(:,1)); ylim([0 10]);
    figure; histogram(cazymedata.M1_barleybran.U{2}(:,2)); ylim([0 10]);
 
    cazymedata.taxa_projections_barleybran.TC1 = cazymedata.M1_barleybran.U{2}(:,1);
    cazymedata.taxa_projections_barleybran.TC2 = cazymedata.M1_barleybran.U{2}(:,2);
   
 
