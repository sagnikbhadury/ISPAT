# ISPAT : Informed Spatially aware Patterns for Multiplexed Immunofluorescence Data
with multi modal integration for Poisson and Gaussian Spatial mixed effects model is possible with this repository

```
  mydata <- readRDS(file = paste0(savepath,"data/",
                                  "Data - Analysis ready Marginally obtained Spatial intensities of Cells on Tissue Slide_", diseasetype,"_",mypatID,".rds"), refhook = NULL)
  rawdata <- mydata %>% pivot_wider(names_from = Cells, values_from = Intensities) #%>% as.matrix()
  colnames(rawdata)
  S <- rawdata[,c("X","Y","area_category")] %>% as.matrix()
  Y <- rawdata %>% select("APC", "CTL", "Epi", "THelper", "Treg") %>% t() %>% as.matrix()
  time<-proc.time()
  saved<-ISPAT(Y, S, ncores = numcores, RefPrior = diag(1,5,5), use_ref = FALSE,
                                Kernel = "Matern",sGLM_method = "MLE", VB_MSFA= TRUE, MSFA_method = "CAVI")
  time<-proc.time()-time
```
