

TTdat <- readRDS(file="../ProcessedData/Combined_tracks_RNAi_2.Rds")
ff <- "2021-07-02_RNAi-3954-3_group1_26"



tt1 <- TTdat[TTdat$id==ff,]
tt1 <- tt1[tt1$likelihood >= 0.6,]



plot(tt1$Second, tt1$x)



