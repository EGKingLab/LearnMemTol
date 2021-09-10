

TTdat <- readRDS(file="../ProcessedData/Combined_tracks_RNAi_2.Rds")
ff <- "2021-06-16_RNAi-50656.25750-1_group1_21"



tt1 <- TTdat[TTdat$id==ff,]
tt1 <- tt1[tt1$likelihood >= 0.6,]



plot(tt1$Second, tt1$x)

ggplot(tt1, aes(Second, x, color=likelihood)) +
  geom_point()



