# Convert Rmd to R script

CompLoc = "office"
#CompLoc = "home"

if(CompLoc=="office")
 {
  knitr::purl(input = "C:\\Users\\weixing\\Dropbox\\LmmSimul\\codes\\Multi-T.Rmd", 
              output = "C:\\Users\\weixing\\Dropbox\\LmmSimul\\codes\\Multi-T.R",
              documentation = 0)

  knitr::purl(input = "C:\\Users\\weixing\\Dropbox\\LmmSimul\\codes\\Multi-Laplace.Rmd", 
              output = "C:\\Users\\weixing\\Dropbox\\LmmSimul\\codes\\Multi-Laplace.R",
              documentation = 0)
  
  knitr::purl(input = "C:\\Users\\weixing\\Dropbox\\LmmSimul\\codes\\Ldegree.Rmd", 
              output = "C:\\Users\\weixing\\Dropbox\\LmmSimul\\codes\\Ldegree.R",
              documentation = 0)
 } else
if(CompLoc=="home")
 {
  knitr::purl(input = "C:\\Users\\weixi\\Dropbox\\LmmSimul\\codes\\Multi-T.Rmd", 
              output = "C:\\Users\\weixi\\Dropbox\\LmmSimul\\codes\\Multi-T.R",
              documentation = 0)
  
  knitr::purl(input = "C:\\Users\\weixi\\Dropbox\\LmmSimul\\codes\\Multi-Laplace.Rmd", 
              output = "C:\\Users\\weixi\\Dropbox\\LmmSimul\\codes\\Multi-Laplace.R",
              documentation = 0)
  
  knitr::purl(input = "C:\\Users\\weixi\\Dropbox\\LmmSimul\\codes\\Ldegree.Rmd", 
              output = "C:\\Users\\weixi\\Dropbox\\LmmSimul\\codes\\Ldegree.R",
              documentation = 0)
 }