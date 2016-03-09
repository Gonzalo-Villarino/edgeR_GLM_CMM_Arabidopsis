#### edgeR_GLM_CMM_Arabidopsis

###If you have more than 3 samples you can make a model. The model knows a bit more about the dataset b/c you provide all the data at once. So you can do more comparisons to be made with statiscal signifiancy, for exampple:

###my.contrasts <- makeContrasts(
        c1 = (PRTYFP_POS-PRTYFP_NEG), # comparison: yfp pos / yfp neg
        c2 = (PRTALL_SORT-PRTNO_SORT),  # comparison all sort / no sort 
        c3 = (PRTALL_SORT-PRTYFP_POS), # comparison all sort / yfp pos
        c4 = (PRTALL_SORT-PRTYFP_NEG), # comparison all sort / yfp neg
        c5 = (PRTNO_SORT-PRTYFP_POS), # comparison no sort / yfp pos
        c6 = (PRTNO_SORT-PRTYFP_NEG), levels=design # comparison no sort / yfp neg
)

You can also ask the difference in genes between c1 and c2.  

what is teh prob. of a CTR gene 
