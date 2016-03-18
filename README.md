#### edgeR_GLM_CMM_Arabidopsis

##### If you have more than 3 samples you can make a model. The model knows a bit more about the dataset b/c you provide all the data at once. So you can do more comparisons to be made with statiscal signifiancy, for exampple, you can also ask the difference in genes between c1 and c2 (c7):

        my.contrasts <- makeContrasts(
        c1 = (PRTYFP_POS-PRTYFP_NEG), # comparison: yfp pos / yfp neg
        c2 = (PRTALL_SORT-PRTNO_SORT),  # comparison all sort / no sort 
        c3 = (PRTALL_SORT-PRTYFP_POS), # comparison all sort / yfp pos
        ...
        ...
        ...
        c7 = (PRTYFP_POS-PRTYFP_NEG)-(PRTALL_SORT-PRTNO_SORT), levels=design # comparison 
        )
        interesting.contrasts <- c("c1", "c2", "c3", "c4", "c5", "c6","c7")

          C7=this comparison looks to see how DEGs behaive differently betwn YFP+/YFP- vs All/NO-sort

##### what is a DEG?
##### What is the probability of a gene in a given condition whose expresiuon value is centered at 10 to be found within the distribution of the same gene in other condition 

