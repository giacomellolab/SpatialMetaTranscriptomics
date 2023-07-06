library('eulerr')

PS <- read.csv2("~/Desktop/Arabidopsis/Figures_data/Supplementary/ExFig1-lableafProfiles/230312_venn_Pseudomonas_C1.csv", sep = ",")
PS <- PS[,-1]
PR1 <- read.csv2("~/Desktop/Arabidopsis/Figures_data/Supplementary/ExFig1-lableafProfiles/230312_venn_PR1_C1.csv", sep = ",")
PR1 <- PR1 [,-1]
Flmax <- read.csv2("~/Desktop/Arabidopsis/Figures_data/Supplementary/ExFig1-lableafProfiles/230312_venn_FLmax_C1.csv", sep = ",")
Flmax <- Flmax[,-1]


# making the intersections
PS_PR1 <- intersect(PS, PR1)
PS_Flmax <- intersect(PS, Flmax)
PR1_Flmax <- intersect(PR1, Flmax)

PS_PR1_Flmax <- intersect(intersect(PS, PR1), Flmax)

PS_only <- setdiff(setdiff(PS, PR1), Flmax)
PR1_only <- setdiff(setdiff(PR1, PS), Flmax)
Flmax_only <- setdiff(setdiff(Flmax, PS), PR1)

PR1_FL_noPS <- setdiff(intersect(Flmax,PR1), PS)

fit1 <- euler(c("PS" = length(PS_only), "PR1" = length(PR1_only), "Flmax" = length(Flmax_only),
                "PS&PR1" = length(PS_PR1), "PS&Flmax" = length(PS_Flmax), "PR1&Flmax" = length(PR1_Flmax),
                "PS&PR1&Flmax" = length(PS_PR1_Flmax)))

#! pdf("~/PhD_stuff/Euller_PS_PR1_Fl.pdf", width = 10, height = 10)
plot(fit1, fills = c("#7fc97f", "#beaed4", "#fdc086"))
#! dev.off()

writeLines(capture.output(sessionInfo()), "~/Desktop/Arabidopsis/Revision_NBT/Rcode/sessionInfo_230706.txt")
