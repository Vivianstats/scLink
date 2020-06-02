lowerTriInd <- function(p){

  i1  <- rep(0,p*(p-1)/2)
  i2  <- rep(0,p*(p-1)/2)
  count <- 1
  for (j in 1:(p-1))
  {
   for (i in (j+1):p){
     i1[count] <- j
     i2[count] <- i
     count <- count+1
   }
  }
 return(list(i1,i2))

}


lowerTriMatInd<- function(x, P){
  i1 <- ceiling(((P-1/2)-sqrt((P-1/2)^2 -2*x)))
  i2 <- P-(i1*P-i1*(i1+1)/2-x)
  return(list(i1,i2))
}