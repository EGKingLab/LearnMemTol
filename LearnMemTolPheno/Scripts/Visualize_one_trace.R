

TTdat <- readRDS(file="../ProcessedData/Combined_tracks_RNAi_2.Rds")
ff <- "2021-03-23_RNAi-51511.51941-1_group1_17"



tt1 <- TTdat[TTdat$id==ff,]
tt1 <- tt1[tt1$likelihood >= 0.6,]



plot(tt1$Second, tt1$x)

ggplot(tt1, aes(Second, x, color=likelihood)) +
  geom_point()
  


