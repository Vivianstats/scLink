########################## Graph structure generation ##########################
graphStructure <- function(nnodes, low.strength, sup.strength, pattern ="hubs", nhubs, 
                           alpha = 2.3, degree.hubs, nOtherEdges, orderNodes=0, 
                           seed = 200, plus = 0, prob=0.05, probSign=0.5)
{
  p           <- nnodes  
  HUBS.NODES  <- 0
  NON.CON.EDG <- 0
  
  if(pattern == "hubs")
  {
      NOD               <- 1:nnodes
      set.seed(seed+19)
      
      # hub nodes
      HUBS.NODES  <- sample(NOD,nhubs)
      
      # connect hub nodes  with other nodes
      HUBS.CON  <- apply(as.matrix(HUBS.NODES),1,function(x){
        r            <- .Random.seed
        .Random.seed <- r + 19
        sample(NOD[-HUBS.NODES], degree.hubs)
      })
      EdgesToChange   <- do.call(rbind,lapply(as.matrix(1:nhubs),function(i) 
        cbind(rep(HUBS.NODES[i], degree.hubs), HUBS.CON[,i])))
      
      # find edges coming from non-hub nodes randomly
      PASS           <- 1:(nnodes*(nnodes-1)/2) 
      PASS           <- sample(PASS)[1:nOtherEdges]
      IS             <- lowerTriMatInd(PASS,nnodes)
      i1             <- IS[[1]]
      i2             <- IS[[2]]
      EdgesToChange2 <- rbind(EdgesToChange,cbind(i1,i2))
      EdgesToChange1 <- EdgesToChange2 + orderNodes
      NON.CON.EDG    <- NOD[!apply(as.matrix(NOD),1,function(i) any(i==EdgesToChange2))] +  orderNodes  
  }  
  if(pattern == "powerLaw")
  {
      # power law distribution (probability of degree k)
      Pk          <- apply(as.matrix(1:30),1, function(k) (k + plus)^(-alpha)/zeta(alpha)   )
      set.seed(seed+19)
      rm1         <- rmultinom(1, nnodes, prob = Pk)
      DEGREE      <- do.call(c,apply(as.matrix(1:length(rm1)),1,function(i) if (rm1[i]>0) rep(i,rm1[i])))
      degreeVert  <- DEGREE[sample(1:length(DEGREE))]
      VV          <- apply(as.matrix(1:nnodes),1,function(i) rep(i,degreeVert[i]))
      if(class(VV)=="list")   Vert.num    <- do.call(c,VV)
      else   Vert.num    <- VV
      
      if(length(Vert.num)>0){
       stub.num    <- 1:length(Vert.num)
       OP          <- stub.num
       Connections <- array(0,dim=c(length(OP)/2,2))
       co          <- 1
       while( length(OP) > 1 ){
        Connections[co,] <- sample(OP)[1:2]
        OP <- OP[-c(which(OP == Connections[co,1]), which(OP == Connections[co,2]) )]
        co <- co + 1
       }
        
       Vert.conn         <- cbind(Vert.num[Connections[,1]],Vert.num[Connections[,2]])
       correctEdg        <- apply(Vert.conn,1,function(x) x[1]!=x[2])
       EdgesToChange2    <- Vert.conn[correctEdg,] 
       EdgesToChange1    <- EdgesToChange2 + orderNodes
      }
      else
       stop("zero vertices in the graph")
      
  }
  if(pattern == "random")
  {
      set.seed(seed+19)
      rb              <- rbinom(nnodes*(nnodes-1)/2, 1, prob = prob)
      AUX             <- rep(0, nnodes*(nnodes-1)/2)
      EdgesToChange2  <- as.matrix(which(rb==1))
      IS              <- lowerTriMatInd(EdgesToChange2,nnodes)
      i1              <- IS[[1]]
      i2              <- IS[[2]]
      EdgesToChange2  <- cbind(i1,i2)
      EdgesToChange1  <- EdgesToChange2 + orderNodes
  
  }
  
  # Simulate precision matrix
  set.seed(seed)
  SignInEdges     <- rbinom(dim(EdgesToChange2)[1],1,probSign)
  SignInEdges     <- ifelse(SignInEdges==0,-1,1)
  set.seed(seed*3)
  IntenInEdges    <- runif(dim(EdgesToChange2)[1], low.strength, sup.strength) * SignInEdges

  PrecMat                <- array(0,dim = c(nnodes, nnodes))
  PrecMat[EdgesToChange2]<- IntenInEdges
  PrecMat                <- PrecMat + t(PrecMat)
  diag(PrecMat)          <- 1

  # Regularization
  PrecMatPD <- PrecMat
  delta     <- 0
  while(max(eigen(PrecMatPD)$val)/min(eigen(PrecMatPD)$val) > nnodes| min(eigen(PrecMatPD)$val) <0)
  {
      delta     <- delta + 0.01
      PrecMatPD <- PrecMatPD + delta * diag(rep(1,nnodes)) 
  }
  
  return(list(edgesToChange = EdgesToChange1, omega = PrecMatPD, hubs = HUBS.NODES + orderNodes, pattern = pattern))
}


