library(RWeka)
library(igraph)
library(Matrix)
library(ggplot2)
library(poweRlaw)
# if run on Mac use this
# dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_74.jdk/Contents/Home/jre/lib/server/libjvm.dylib')

# harmonic centrality as it is defined in 
# Boldi and Vigna (2013)
harmonic <- function(g){
  
  tmp <- 1/shortest.paths(g,mode ="out");
  
  # computes sums by column, removing infinite values
  myfun <- function(x) {tmp1 <- x[which(x!=Inf)]; return(sum(tmp1));}
  
  # returns sum of reciprocal of distances by column
  return(apply(tmp,2, myfun));
}

# delivers non-discretized values of clues and elapsed times
# that is, the in-degree and the out-degree 
# there is a repetition of these variables since in a previous version
# these variables slots were used to compute other clues
get_raw_clues_2 <- function(g){
  
  # matrix to return
  M <- matrix(nrow=vcount(g),ncol=0);
  
  # get in-degree
  tmp <- proc.time();
  in.deg<-degree(g,mode="in");
  ind.t <- proc.time() - tmp;
  M <- cbind(M, in.deg);
  
  tmp <- proc.time();
  na <- in.deg
  na.t <- proc.time() - tmp;
  M <- cbind(M,na);
  
  tmp <- proc.time();
  bet.pos <- in.deg
  bet.pos.t <- proc.time() - tmp;
  M <- cbind(M, bet.pos);
  
  # get out-degree  
  tmp <- proc.time();
  out.deg<-degree(g,mode="out");
  out.t <- proc.time() - tmp;
  M <- cbind(M, out.deg);
  
  tmp <- proc.time();
  nc <- in.deg
  nc.t <- proc.time() - tmp;
  M <- cbind(M, nc);
  
  tmp <- proc.time();
  bet.neg <- in.deg
  bet.neg.t <- proc.time() - tmp;
  M <- cbind(M, bet.neg);
  
  return(list(M, list(ind.t,
                      na.t,
                      bet.pos.t,
                      out.t,
                      nc.t,
                      bet.neg.t)));
}

# computes distinct norms of the estimation error
# made by the distinct methods (L, NN, QC, T)
# n_iter is the number of networks to produce 
# pv is the vector of proportions (for the quantile degree values)
# prop_oth is the proportion of the dataset used as training set
get_error_by_method_3 <- function(n_iter,pv,prop_oth){
  
  # ncol is the number of methods to compare
  resultados.L1 <- matrix(nrow = n_iter, ncol = 8);
  resultados.L2 <- matrix(nrow = n_iter, ncol = 8);
  resultados.Linf <- matrix(nrow = n_iter, ncol = 8);
  resultados.corr <- matrix(nrow = n_iter, ncol = 8);
  resultados.signif <- matrix(nrow = n_iter, ncol = 8);
  resultados.rmsd <- matrix(nrow = n_iter, ncol = 8);
  resultados.mae <- matrix(nrow = n_iter, ncol = 8);
  
  # tiempos stores the elapsed time of each method
  # 35 = 7 * 5, 5 is the length of each vector returned by proc.time()
  tiempos <- matrix(nrow=n_iter,ncol=35);
  
  for(i in 1:n_iter){
    
    gr <- sample_pa(10000,power=1);
    PT <- get_raw_clues_2(gr);
    results <- get_comparison_3(gr,pv,c(1,0,0,0,0,0,0,0),PT,1,
                                 training_size=prop_oth,tom_prop = prop_oth);
    
    # store results
    
    resultados.L1[i,] <- results[[12]][,1];
    resultados.L2[i,] <- results[[12]][,2];
    resultados.Linf[i,] <- results[[12]][,3];
    resultados.corr[i,] <- results[[12]][,4];
    resultados.signif[i,] <- results[[12]][,5];
    resultados.rmsd[i,] <- results[[12]][,6];
    resultados.mae[i,] <- results[[12]][,7];
    
    # linear, Tree, NN, partI, clues, averages, QC
    tiempos[i,] <- cbind(results[[15]],
                         results[[18]],
                         results[[20]],
                         results[[13]],
                         results[[6]],
                         results[[9]][[2]],
                         results[[14]]);
    
  }
  
  return(list(resultados.L1, 
              resultados.L2, 
              resultados.Linf, 
              resultados.corr, 
              resultados.signif, 
              resultados.rmsd, 
              resultados.mae,
              tiempos));
} 

# receives a sample of values, and computes
# the MLE estimator of the alpha exponent of the distribution
# x_min is the minimum value of the distribution
# if x_min is provided, values smaller than x_min are excluded
get_alpha <- function(values,x_min=0){
  
  if(x_min==0){# x_min was not provided
    x_min = min(values);
    valid_pos <- rep(1,length(values));
  }
  else{# exclude smaller values than x_min
    valid_pos <- which(values >= x_min);
  }
  
  aux <- values[valid_pos];
  n <- length(aux); 
  
  return(1 + n/sum(log(aux/x_min)) );
}

# compute quantile values assuming a power-law
# receives alpha, x_min and vect_prob (vector of quantile probabilities)
get_quantile_6 <- function(alpha,x_min,vect_prob){
  
  return(list(vect_prob,
              x_min * (1 - vect_prob)^ (1/(1 - alpha)) 
  ));
  
}

# delivers discretized clues
# and elapsed computation times
# PT is the output of get_raw_clues()
# c_indir is TRUE iff compute clue quantile values as pre-image by CDF of centrality quantile
# cent_nam, ng are printed in the warning message
# several clues are deprecated
get_final_clues_3 <- function(PT,probs=c(0.1,0.5,0.8),pr_Beta=0.4,
                              instruc=c(1,1,1,1,1,1,1,1),cent,
                              c_indir=TRUE,cent_name="unknown",ng=1){ 
  
  # matrix to return
  M <- matrix(nrow=nrow(PT[[1]]),ncol=0);
  
  total_time <- c(0,0,0,0,0);
  
  # names of the clues to compute
  # these names are deprecated since correspond to variables
  # that were computed in a previous version
  nombres_pistas <- c("IND","REA","BEP","OUT","COR","BEN");
  
  in.deg <- c();
  out.deg <- c();
  
  # look for the indexes in cent, of the exact quantile values given by probs
  if (c_indir) r.cent <- match(quantile(cent,probs,names=FALSE,type=1),
                               sort(cent));
  
  # if there is more than 1 clue selected, compute the product of them
  pos_pistas <- which(instruc == 1);
  if(length(pos_pistas) > 1){
    
    # conversion of clue indexes (deprecated)
    instruc_to_col <- c(1,2,3,3,4,5,6,6);
    cols <- unique(instruc_to_col[pos_pistas]);
    
    tmp <- proc.time();
    
    comb <- PT[[1]][,cols[1]];
    for(i in 2:length(cols)){
      comb <- comb * PT[[1]][,cols[i]];
    }
    
    real_probs <- probs;
    
    # get sample quantile values according to required probabilities
    if (c_indir) {comb.c <- sort(comb)[r.cent]}
    else {comb.c<-quantile(comb,probs,names=FALSE,type=1);}
    
    comb.c.distintos <- unique(comb.c);
    
    # check for repeated quantiles
    try(if(length(comb.c)!=length(comb.c.distintos)){
      
      warning(paste("Quantile vector has repeated elements! With clue comb, \n",
                    "Probabilities =",toString(probs),", \n",
                    "Q. vector =",toString(comb.c),", \n",
                    cent_name," centrality, \n",
                    "n.graph = ", ng));
      
      real_probs <- probs[match(comb.c.distintos,comb.c)];
      comb.c <- comb.c.distintos;
    });
    
    # compute clues, that combination of predictors is greater than respective quantile
    comb.p <- sapply(comb.c, function(x){comb > x});
    
    # sum the total elapsed time, including the computation of non-discretized clues 
    aux <- PT[[2]][[ cols[1] ]];
    for(i in 2:length(cols)){
      aux <- aux + PT[[2]][[ cols[i] ]];
    }
    
    total_time <- total_time + proc.time() - tmp + aux;
    
    # merge clue names
    colnames(comb.p) <- sapply(real_probs, function(x){return(
      paste(paste(nombres_pistas[cols],collapse="_"),":",toString(x),sep=""))});
    
    # M contains all clues
    M <- cbind(M, comb.p);
    
  }
  else{# only 1 clue required
    
    if (instruc[1]==1){# in-degree
      
      in.deg<-PT[[1]][,1];
      
      real_probs <- probs;
      
      tmp <- proc.time();
      
      # get sample quantile values according to required probabilities
      if (c_indir) {in.deg.c <- sort(in.deg)[r.cent]}
      else {in.deg.c<-quantile(in.deg,probs,names=FALSE,type=1);}
      
      in.deg.c.distintos <- unique(in.deg.c);
      
      # check for repeated quantiles
      try(if(length(in.deg.c)!=length(in.deg.c.distintos)){
        
        warning(paste("Quantile vector has repeated elements! With clue IND, \n",
                      "Probabilities =",toString(probs),", \n",
                      "Q. vector =",toString(in.deg.c),", \n",
                      cent_name," centrality, \n",
                      "n.graph = ", ng));
        
        real_probs <- probs[match(in.deg.c.distintos,in.deg.c)];
        in.deg.c <- in.deg.c.distintos;
      });
      
      # compute clues, that in-degree is greater than respective quantile
      in.deg.p <- sapply(in.deg.c, function(x){in.deg > x});
      
      # sum the total elapsed time, including the computation of non-discretized clues 
      total_time <- total_time + proc.time() - tmp + PT[[2]][[1]];
      
      colnames(in.deg.p) <- sapply(real_probs, function(x){return(paste("IND:",toString(x)))})
      
      # M contains all clues
      M <- cbind(M, in.deg.p);
      
    }
    
    if(instruc[2]==1){# number of reachable nodes
      
      na <- PT[[1]][,2]
      
      real_probs <- probs;
      
      tmp <- proc.time();
      
      # get sample quantile values according to required probabilities
      if (c_indir) {na.c <- sort(na)[r.cent]}
      else {na.c<-quantile(na,probs,names=FALSE,type=1);}
      
      na.c.distintos <- unique(na.c);
      
      # check for repeated quantiles
      try(if(length(na.c)!=length(na.c.distintos)){
        
        warning(paste("Quantile vector has repeated elements! With clue REA, \n",
                      "Probabilities =",toString(probs),", \n",
                      "Q. vector =",toString(na.c),", \n",
                      cent_name," centrality, \n",
                      "n.graph = ", ng));
        
        real_probs <- probs[match(na.c.distintos,na.c)];
        na.c <- na.c.distintos;
      });
      
      # compute clues, that number of reachable nodes is greater than respective quantile
      na.p <- sapply(na.c, function(x){na > x});
      
      # sum the total elapsed time, including the computation of non-discretized clues 
      total_time <- total_time + proc.time() - tmp + PT[[2]][[2]];
      
      colnames(na.p) <- sapply(real_probs, function(x){return(paste("REA:",toString(x)))})
      
      M <- cbind(M,na.p);
    }
    
    if(instruc[3]==1){ # slot for a previous variable non-valid currently (deprecated)
      
      in.deg <- PT[[1]][,1];
      
      tmp <- proc.time();
      
      be.p <- sapply(probs, function(x) {Beta_positivo(g,
                                                       in.deg <= quantile(in.deg,pr_Beta,names=FALSE),
                                                       thre = x);}
      );
      
      total_time <- total_time + proc.time() - tmp + PT[[2]][[1]];
      
      colnames(be.p) <- sapply(probs, function(x){return(paste("BP-",toString(x)))})
      
      M <- cbind(M, be.p);
    }
    
    if(instruc[4]==1){# slot for a previous variable non-valid currently (deprecated)
      
      bet.pos <- PT[[1]][,3];
      
      real_probs <- probs;
      
      tmp <- proc.time();
      
      if (c_indir) {bet.pos.c <- sort(bet.pos)[r.cent]}
      else {bet.pos.c<-quantile(bet.pos,probs,names=FALSE,type=1);}
      
      bet.pos.c.distintos <- unique(bet.pos.c);
      
      try(if(length(bet.pos.c)!=length(bet.pos.c.distintos)){
        
        warning(paste("Quantile vector has repeated elements! With clue BEP, \n",
                      "Probabilities =",toString(probs),", \n",
                      "Q. vector =",toString(bet.pos.c),", \n",
                      cent_name," centrality, \n",
                      "n.graph = ", ng));
        
        real_probs <- probs[match(bet.pos.c.distintos,bet.pos.c)];
        bet.pos.c <- bet.pos.c.distintos;
      });

      bet.pos.p <- sapply(bet.pos.c, function(x){bet.pos > x});
      
      total_time <- total_time + proc.time() - tmp + PT[[2]][[3]] + PT[[2]][[1]];
      
      colnames(bet.pos.p) <- sapply(real_probs, function(x){return(paste("BEP:",toString(x)))})
      
      M <- cbind(M, bet.pos.p);
    }
    
    if (instruc[5]==1){# out-degree
      
      out.deg<-PT[[1]][,4];
      
      real_probs <- probs;
      
      tmp <- proc.time();
      
      # get sample quantile values according to required probabilities
      if (c_indir) {out.deg.c <- sort(out.deg)[r.cent]}
      else {out.deg.c<-quantile(out.deg,probs,names=FALSE,type=1);}
      
      out.deg.c.distintos <- unique(out.deg.c);
      
      # check for repeated quantiles
      try(if(length(out.deg.c)!=length(out.deg.c.distintos)){
        
        warning(paste("Quantile vector has repeated elements! With clue OUT, \n",
                      "Probabilities =",toString(probs),", \n",
                      "Q. vector =",toString(out.deg.c),", \n",
                      cent_name," centrality, \n",
                      "n.graph = ", ng));
        
        real_probs <- probs[match(out.deg.c.distintos,out.deg.c)];
        out.deg.c <- out.deg.c.distintos;
      });
      
      # compute clues, that out-degree is greater than respective quantile
      out.deg.p <- sapply(out.deg.c, function(x){out.deg > x});
      
      # sum the total elapsed time, including the computation of non-discretized clues 
      total_time <- total_time + proc.time() - tmp + PT[[2]][[4]];
      
      colnames(out.deg.p) <- sapply(real_probs, function(x){return(paste("OUT:",toString(x)))})
      
      # M contains all clues
      M <- cbind(M, out.deg.p);
    }
    
    if(instruc[6]==1){# number of coreachable nodes
      
      nc <- PT[[1]][,5];
      
      real_probs <- probs;
      
      tmp <- proc.time();
      
      # get sample quantile values according to required probabilities
      if (c_indir) {nc.c <- sort(nc)[r.cent]}
      else {nc.c<-quantile(nc,probs,names=FALSE,type=1);}
      
      nc.c.distintos <- unique(nc.c);
      
      # check for repeated quantiles
      try(if(length(nc.c)!=length(nc.c.distintos)){
        
        warning(paste("Quantile vector has repeated elements! With clue COR, \n",
                      "Probabilities =",toString(probs),", \n",
                      "Q. vector =",toString(nc.c),", \n",
                      cent_name," centrality, \n",
                      "n.graph = ", ng));
        
        real_probs <- probs[match(nc.c.distintos,nc.c)];
        nc.c <- nc.c.distintos;
      });
      
      # compute clues, that number of coreachable is greater than respective quantile
      nc.p <- sapply(nc.c, function(x){nc > x});
      
      # sum the total elapsed time, including the computation of non-discretized clues 
      total_time <- total_time + proc.time() - tmp + PT[[2]][[5]];
      
      colnames(nc.p) <- sapply(real_probs, function(x){return(paste("COR:",toString(x)))})
      
      # M contains all clues
      M <- cbind(M, nc.p);
    }
    
    if(instruc[7]==1){ # slot for a previous variable non-valid currently (deprecated)
      
      out.deg <- PT[[1]][,4];
      
      tmp <- proc.time();
      
      be.n.p <- sapply(probs, function(x) {Beta_negativo(g,
                                                         out.deg <= quantile(out.deg,pr_Beta,names=FALSE),
                                                         thre=x);}
      );
      
      total_time <- total_time + proc.time() - tmp + PT[[2]][[4]];
      
      colnames(be.n.p) <- sapply(probs, function(x){return(paste("BN:",toString(x)))})
      
      M <- cbind(M, be.n.p);
    }
    
    if(instruc[8]==1){# slot for a previous variable non-valid currently (deprecated)
      
      bet.neg <- PT[[1]][,6];
      
      real_probs <- probs;
      
      tmp <- proc.time();
      
      if (c_indir) {bet.neg.c <- sort(bet.neg)[r.cent]}
      else {bet.neg.c<-quantile(bet.neg,probs,names=FALSE,type=1);}
      
      bet.neg.c.distintos <- unique(bet.neg.c);
      
      try(if(length(bet.neg.c)!=length(bet.neg.c.distintos)){
        
        warning(paste("Quantile vector has repeated elements! With clue BEN, \n",
                      "Probabilities =",toString(probs),", \n",
                      "Q. vector =",toString(bet.neg.c),", \n",
                      cent_name," centrality, \n",
                      "n.graph = ", ng));
        
        real_probs <- probs[match(bet.neg.c.distintos,bet.neg.c)];
        bet.neg.c <- bet.neg.c.distintos;
      });
      
      bet.neg.p <- sapply(bet.neg.c, function(x){bet.neg > x});
      
      total_time <- total_time + proc.time() - tmp + PT[[2]][[6]] + PT[[2]][[4]];
      
      colnames(bet.neg.p) <- sapply(real_probs, function(x){return(paste("BEN:",toString(x)))})
      
      M <- cbind(M, bet.neg.p);
      
    }
  } # else end
  
  return(list(M, total_time));
}

# compute the representative median associated to each interval 
# returns an object containing the vector of sorted (ascending) medians 
# and the sorting to traverse the clue columns
# c_quant is the vector quantile centrality values
# alfa is the exponent of the power law
# x_min the minimum value of this distribution
get_Prom_Ord_5 <- function(c_quant,alfa,x_min,cl_instr,pistas){
  
  # clue names (deprecated)
  clue_names <- c("IND","REA","BP","BEP","OUT","COR","BN","BEN");
  # clues used to build the clue matrix
  current_clues <- clue_names[which(cl_instr==1)];
  # name for each column
  clue_names <- unlist(strsplit(colnames(pistas), ":", fixed = TRUE))[seq(from=1,to=dim(pistas)[2]*2,by=2)];
  # object that stores the median of each clue
  promedios.pistas <- rep(0,dim(pistas)[2] + length(current_clues));
  
  # column names
  name_pr <- rep(0,dim(pistas)[2] + length(current_clues));
  
  total_time <- c(0,0,0,0,0);
  
  # multiplicative constant
  C_mult <- 2^(1/(alfa - 1));
  
  for(i in seq_along(current_clues)){
    
    clue_indexes <- grep(current_clues[i],clue_names)
    
    time_stamp_1 <- proc.time();
    
    # median of the first set
    # assumes that probability quantiles in c_quant comes in ascending order
    promedios.pistas[clue_indexes[1]] <- C_mult * (c_quant[1]^(1-alfa) + x_min^(1-alfa))^(1/(1 - alfa));
    
    total_time = total_time + proc.time() - time_stamp_1;
    
    # set column name
    name_pr[clue_indexes[1]] <- colnames(pistas)[clue_indexes[1]];
    
    # if there are more columns based on the same clue (and distinct quantile)
    if(length(clue_indexes)>1){
      for (j  in 2:length(clue_indexes)){
        
        time_stamp <- proc.time();
        
        promedios.pistas[clue_indexes[j]] <- C_mult * (c_quant[j]^(1-alfa) + c_quant[j-1]^(1-alfa))^(1/(1 - alfa));
        
        total_time = total_time + proc.time() - time_stamp;
        
        # set column name
        name_pr[clue_indexes[j]] <- colnames(pistas)[clue_indexes[j]];
        
      }
    }
    
    # median of the last set
    
    time_stamp_3 <- proc.time();
    
    promedios.pistas[dim(pistas)[2] + i] <- C_mult * c_quant[length(clue_indexes)];
    
    total_time = total_time + proc.time() - time_stamp_3;
    
    # set column name
    name_pr[dim(pistas)[2] + i] <- paste(colnames(pistas)[clue_indexes[length(clue_indexes)]],"+1");
    
  }
  
  # assign names to set medians
  names(promedios.pistas) <- name_pr;
  
  time_stamp_2 <- proc.time();
  
  # now promedios.pistas has 2 fields, $x, medians in ascending order, 
  # and $ix, the order in which to traverse columns
  promedios.pistas <- sort(promedios.pistas,index.return = TRUE);
  
  total_time = total_time + proc.time() - time_stamp_2;
  
  return(list(promedios.pistas,total_time));
}

# QuickCent method, returns clue index used to get the estimate
QuickCent <- function(prom_ord,i,pistas){
  
  # clues of the i-th object  
  pistas_i <- pistas[i,]; 
  
  # get indexes of false clues
  index_false <- which(!pistas_i)
  
  # get the least index from prom_ord$ix
  # where there is a match with index_false
  # first case is when there are no false clues
  if (length(index_false)==0)
    return(length(pistas_i) + 1);
  
  return(min(match(index_false,prom_ord$ix),na.rm = TRUE));
  
}

# Compares performance of QC to other ML methods 
# returns estimation error norms 
# receives graph g, vector of quantile proportions vect_prob 
# cl_instr is the set of clues
# MC is the output of get_raw_clues_2()
# builds clue matrices with get_final_clues_3()
# computes averages with get_Prom_ord_5()
# n_graph es el número de grafo que estamos procesando
get_comparison_3 <- function(g,vect_prob,
                              cl_instr = c(1,0,0,0,0,0,0,0),
                              MC,n_graph,
                              x_min=1,
                              training_size=0.2,
                              tom_prop=0.3
){
  
  # compute harmonic centrality (part I)
  c_time_1 <- proc.time();
  harm.cent <- harmonic(g);
  c_time_2 <- proc.time();
  cent_comput_time <- c_time_2 - c_time_1; 
  
  time_partI <- proc.time();
  
  # estimate alpha by sampling from the total of nodes set
  total_nodes <- length(harm.cent);
  alfa <- get_alpha(harm.cent[sample(1:total_nodes,
                                     trunc(total_nodes*training_size))],x_min);
  
  harm.cuant <- get_quantile_6(alfa,x_min,vect_prob);
  
  time_partI <- proc.time() - time_partI;
  
  # get discretized clues
  results.harm <- get_final_clues_3(MC,probs=unique(vect_prob),instruc= cl_instr,
                                    cent =harm.cent,c_indir=FALSE,cent_name="Harmonic",ng=n_graph)
  
  #pistas.harm has, at each column, the values of one clue
  pistas.harm <- results.harm[[1]];
  t_pistas.harm <- results.harm[[2]];
  
  # compute average of each clue, get clues sorting according to their averages
  
  p.ord.harm <- get_Prom_Ord_5(harm.cuant[[2]],alfa,x_min,cl_instr,pistas.harm);
  
  # get QuickCent estimate (part II)
  time_partII <- proc.time();
  
  est.harm.index <- sapply(c(1:total_nodes),function(i) {QuickCent(p.ord.harm[[1]],i,pistas.harm)})
  est.harm <- p.ord.harm[[1]]$x[est.harm.index];
  
  time_partII <- proc.time() - time_partII;
  
  
  # compute norm of error made by centrality estimate
  
  normas.qcent <- compute_error_norms(est.harm,harm.cent);
  
  # build dataset matrix to apply other regression methods
  deg=MC[[1]][,1];
  df.har <- data.frame(harmonic=harm.cent,
                       pistas.harm,
                       deg);
  # rows_to_train is the sampling used as training set
  rows_to_train <- sample(1:total_nodes,trunc(total_nodes*tom_prop));
  
  # compute linear regression
  
  time_lin.1 <- proc.time();
  
  lin.1 <- lm(harmonic ~  I(deg) , data = df.har[rows_to_train,]);
  pred.lin.1 <- predict(lin.1,newdata = df.har);
  
  time_lin.1 <- proc.time() - time_lin.1;
  
  
  time_lin.2 <- proc.time();
  
  lin.2 <- lm(as.formula(paste('harmonic~',paste(colnames(df.har)[2:(ncol(df.har)-1)],
                                                 collapse='+'))), data = df.har[rows_to_train,]);
  pred.lin.2 <- predict(lin.2,newdata = df.har);
  
  time_lin.2 <- proc.time() - time_lin.2;
  
  
  time_lin.3 <- proc.time();
  
  lin.3 <- lm(as.formula(paste('I(log(harmonic + 1))~',
                               paste(colnames(df.har)[2:(ncol(df.har)-1)],collapse='+'),
                               '+ I(log(deg + 1))')),
              data = df.har[rows_to_train,]);
  pred.lin.3 <- exp(predict(lin.3,newdata = df.har)) - 1;
  
  time_lin.3 <- proc.time() - time_lin.3;
  
  normas.lin.1 <- compute_error_norms(harm.cent,pred.lin.1);
  normas.lin.2 <- compute_error_norms(harm.cent,pred.lin.2);
  normas.lin.3 <- compute_error_norms(harm.cent,pred.lin.3);
  
  # compute Tree and NN classifiers
  
  f.REPTree <- make_Weka_classifier("weka/classifiers/trees/REPTree", c("bar","Weka_tree"))
  #f.SVM <- make_Weka_classifier("weka/classifiers/functions/SMOreg", c("bar","Weka_tree"))
  f.NN <- make_Weka_classifier("weka/classifiers/functions/MultilayerPerceptron", c("bar","Weka_tree"))
  
  time_tree.1 <- proc.time();
  
  tree.1 <- f.REPTree(harmonic ~  I(deg),data = df.har[rows_to_train,]);
  pred.tree.1 <- predict(tree.1,newdata = df.har);
  
  time_tree.1 <- proc.time() - time_tree.1;
  
  
  time_tree.2 <- proc.time();
  
  tree.2 <- f.REPTree(as.formula(paste('harmonic~',paste(colnames(df.har)[2:ncol(df.har)],
                                                         collapse='+'),'+ I(deg)')),data = df.har[rows_to_train,]);
  pred.tree.2 <- predict(tree.2,newdata = df.har);
  
  time_tree.2 <- proc.time() - time_tree.2;
  
  
  normas.tree.1 <- compute_error_norms(harm.cent,pred.tree.1);
  normas.tree.2 <- compute_error_norms(harm.cent,pred.tree.2);
  
  
  time_NN.1 <- proc.time();
  
  NN.1 <- f.NN(harmonic ~  I(deg),data = df.har[rows_to_train,]);
  pred.NN.1 <- predict(NN.1,newdata = df.har);
  
  time_NN.1 <- proc.time() - time_NN.1;
  
  
  time_NN.2 <- proc.time();
  
  NN.2 <- f.NN(as.formula(paste('harmonic~',paste(colnames(df.har)[2:ncol(df.har)],
                                                  collapse='+'),'+ I(deg)')),data = df.har[rows_to_train,]);
  pred.NN.2 <- predict(NN.2,newdata = df.har);
  
  time_NN.2 <- proc.time() - time_NN.2;
  
  
  normas.NN.1 <- compute_error_norms(harm.cent,pred.NN.1);
  normas.NN.2 <- compute_error_norms(harm.cent,pred.NN.2);
  
  normas.harm <- rbind(normas.qcent,
                       normas.lin.1,
                       normas.lin.2,
                       normas.lin.3,
                       normas.tree.1,
                       normas.tree.2,
                       normas.NN.1,
                       normas.NN.2);
  
  data <- list(harm.cent,
               alfa,
               cent_comput_time,
               harm.cuant,
               pistas.harm,
               t_pistas.harm,
               vect_prob,
               cl_instr,
               p.ord.harm,
               est.harm.index,
               # 11
               est.harm,
               normas.harm,
               time_partI,
               time_partII,
               time_lin.1,
               time_lin.2,
               time_lin.3,
               time_tree.1,
               time_tree.2,
               #20
               time_NN.1,
               time_NN.2
  );
  
  return(data);
}

# computes the error made when trained with 
# distinct proportions of the dataset
# n_iter is the number of networks to test
# pv is the vector of quantile probabilities
get_error_by_prop <- function(n_iter,pv){
  
  resultados.L1 <- matrix(nrow = n_iter, ncol = 10);
  resultados.L2 <- matrix(nrow = n_iter, ncol = 10);
  resultados.Linf <- matrix(nrow = n_iter, ncol = 10);
  resultados.corr <- matrix(nrow = n_iter, ncol = 10);
  resultados.signif <- matrix(nrow = n_iter, ncol = 10);
  resultados.rmsd <- matrix(nrow = n_iter, ncol = 10);# 6 is RMSD
  resultados.mae <- matrix(nrow = n_iter, ncol = 10);# 7 is MAE (mean absolute error)
  
  for(i in 1:n_iter){
    
    gr <- sample_pa(10000,power=1);
    PT <- get_raw_clues_2(gr);
    results <- get_simul_6(gr,vect_prob=pv,
                           cl_instr = c(1,0,0,0,0,0,0,0),
                           PT,i,x_min=1
    )
    
    # guardamos los resultados
    resultados.L1[i,] <- results[[58]][,1];
    resultados.L2[i,] <- results[[58]][,2];
    resultados.Linf[i,] <- results[[58]][,3];
    resultados.corr[i,] <- results[[58]][,4];
    resultados.signif[i,] <- results[[58]][,5];
    resultados.rmsd[i,] <- results[[58]][,6];
    resultados.mae[i,] <- results[[58]][,7];
    
  }
  
  return(list(resultados.L1, 
              resultados.L2, 
              resultados.Linf, 
              resultados.corr, 
              resultados.signif, 
              resultados.rmsd, 
              resultados.mae));
} 

# computes norm of estimation error
# receives a graph g, vector of quantile proportions vect_prob
# cl_instr is the set of clues
# MC is the output of get_raw_clues_2()
# builds matrices of clues with get_final_clues_3()
# uses get_Prom_ord_5() to compute representative medians
# n_graph is the number of graph to be processed
get_simul_6 <- function(g,vect_prob,
                        cl_instr = c(1,0,0,0,0,0,0,0),
                        MC,n_graph,
                        x_min=1,
                        training_sizes=seq(1,0.1,-0.1)
){
  
  # compute harmonic centrality
  c_time_1 <- proc.time();
  harm.cent <- harmonic(g);
  c_time_2 <- proc.time();
  cent_comput_time <- c_time_2 - c_time_1; 
  
  # estimation of alpha exponent by sampling from the total of nodes set
  total_nodes <- length(harm.cent);
  alfa <- sapply(training_sizes,
                 function(x){
                   get_alpha(harm.cent[sample(1:total_nodes,
                                              trunc(total_nodes*x))],x_min)});
  
  harm.cuant.1 <- get_quantile_6(alfa[1],x_min,vect_prob);
  harm.cuant.2 <- get_quantile_6(alfa[2],x_min,vect_prob);
  harm.cuant.3 <- get_quantile_6(alfa[3],x_min,vect_prob);
  harm.cuant.4 <- get_quantile_6(alfa[4],x_min,vect_prob);
  harm.cuant.5 <- get_quantile_6(alfa[5],x_min,vect_prob);
  harm.cuant.6 <- get_quantile_6(alfa[6],x_min,vect_prob);
  harm.cuant.7 <- get_quantile_6(alfa[7],x_min,vect_prob);
  harm.cuant.8 <- get_quantile_6(alfa[8],x_min,vect_prob);
  harm.cuant.9 <- get_quantile_6(alfa[9],x_min,vect_prob);
  harm.cuant.10 <- get_quantile_6(alfa[10],x_min,vect_prob);

  # get dicretized clues
  results.harm <- get_final_clues_3(MC,probs=unique(vect_prob),instruc= cl_instr,
                                    cent =harm.cent,c_indir=FALSE,cent_name="Harmonic",ng=n_graph)
  
  # pistas.harm has, at each column, the values of one clue
  pistas.harm <- results.harm[[1]];
  t_pistas.harm <- results.harm[[2]];
  
  # compute sorted set of medians associated to each clue
  # for each sampling size of the dataset 
  
  # 1
  p.ord.harm.1 <- get_Prom_Ord_5(harm.cuant.1[[2]],alfa[1],x_min,cl_instr,pistas.harm);
  # 0.9
  p.ord.harm.2 <- get_Prom_Ord_5(harm.cuant.2[[2]],alfa[2],x_min,cl_instr,pistas.harm);
  # 0.8
  p.ord.harm.3 <- get_Prom_Ord_5(harm.cuant.3[[2]],alfa[3],x_min,cl_instr,pistas.harm);
  # 0.7
  p.ord.harm.4 <- get_Prom_Ord_5(harm.cuant.4[[2]],alfa[4],x_min,cl_instr,pistas.harm);
  # 0.6
  p.ord.harm.5 <- get_Prom_Ord_5(harm.cuant.5[[2]],alfa[5],x_min,cl_instr,pistas.harm);
  # 0.5
  p.ord.harm.6 <- get_Prom_Ord_5(harm.cuant.6[[2]],alfa[6],x_min,cl_instr,pistas.harm);
  # 0.4
  p.ord.harm.7 <- get_Prom_Ord_5(harm.cuant.7[[2]],alfa[7],x_min,cl_instr,pistas.harm);
  # 0.3
  p.ord.harm.8 <- get_Prom_Ord_5(harm.cuant.8[[2]],alfa[8],x_min,cl_instr,pistas.harm);
  # 0.2
  p.ord.harm.9 <- get_Prom_Ord_5(harm.cuant.9[[2]],alfa[9],x_min,cl_instr,pistas.harm);
  # 0.1
  p.ord.harm.10 <- get_Prom_Ord_5(harm.cuant.10[[2]],alfa[10],x_min,cl_instr,pistas.harm);
  
  # get QuickCent estimate
  c_time_4 <- proc.time();
  
  est.harm.index.1 <- sapply(c(1:total_nodes),function(i) {QuickCent(p.ord.harm.1[[1]],i,pistas.harm)})
  est.harm.1 <- p.ord.harm.1[[1]]$x[est.harm.index.1];
  
  c_time_5 <- proc.time();
  
  est.harm.index.2 <- sapply(c(1:total_nodes),function(i) {QuickCent(p.ord.harm.2[[1]],i,pistas.harm)})
  est.harm.2 <- p.ord.harm.2[[1]]$x[est.harm.index.2];
  
  c_time_6 <- proc.time();
  
  est.harm.index.3 <- sapply(c(1:total_nodes),function(i) {QuickCent(p.ord.harm.3[[1]],i,pistas.harm)})
  est.harm.3 <- p.ord.harm.3[[1]]$x[est.harm.index.3];
  
  c_time_7 <- proc.time();
  
  est.harm.index.4 <- sapply(c(1:total_nodes),function(i) {QuickCent(p.ord.harm.4[[1]],i,pistas.harm)})
  est.harm.4 <- p.ord.harm.4[[1]]$x[est.harm.index.4];
  
  c_time_8 <- proc.time();
  
  est.harm.index.5 <- sapply(c(1:total_nodes),function(i) {QuickCent(p.ord.harm.5[[1]],i,pistas.harm)})
  est.harm.5 <- p.ord.harm.5[[1]]$x[est.harm.index.5];
  
  c_time_9 <- proc.time();
  
  est.harm.index.6 <- sapply(c(1:total_nodes),function(i) {QuickCent(p.ord.harm.6[[1]],i,pistas.harm)})
  est.harm.6 <- p.ord.harm.6[[1]]$x[est.harm.index.6];
  
  c_time_10 <- proc.time();
  
  est.harm.index.7 <- sapply(c(1:total_nodes),function(i) {QuickCent(p.ord.harm.7[[1]],i,pistas.harm)})
  est.harm.7 <- p.ord.harm.7[[1]]$x[est.harm.index.7];
  
  c_time_11 <- proc.time();
  
  est.harm.index.8 <- sapply(c(1:total_nodes),function(i) {QuickCent(p.ord.harm.8[[1]],i,pistas.harm)})
  est.harm.8 <- p.ord.harm.8[[1]]$x[est.harm.index.8];
  
  c_time_12 <- proc.time();
  
  est.harm.index.9 <- sapply(c(1:total_nodes),function(i) {QuickCent(p.ord.harm.9[[1]],i,pistas.harm)})
  est.harm.9 <- p.ord.harm.9[[1]]$x[est.harm.index.9];
  
  c_time_13 <- proc.time();
  
  est.harm.index.10 <- sapply(c(1:total_nodes),function(i) {QuickCent(p.ord.harm.10[[1]],i,pistas.harm)})
  est.harm.10 <- p.ord.harm.10[[1]]$x[est.harm.index.10];
  
  c_time_14 <- proc.time();
  
  # compute estimation error
  err.harm.1 <- est.harm.1 - harm.cent;
  err.harm.2 <- est.harm.2 - harm.cent;
  err.harm.3 <- est.harm.3 - harm.cent;
  err.harm.4 <- est.harm.4 - harm.cent;
  err.harm.5 <- est.harm.5 - harm.cent;
  err.harm.6 <- est.harm.6 - harm.cent;
  err.harm.7 <- est.harm.7 - harm.cent;
  err.harm.8 <- est.harm.8 - harm.cent;
  err.harm.9 <- est.harm.9 - harm.cent;
  err.harm.10 <- est.harm.10 - harm.cent;
  
  # variables to store error norms
  normas.harm.1 <- rep(0,7);
  normas.harm.2 <- rep(0,7);
  normas.harm.3 <- rep(0,7);
  normas.harm.4 <- rep(0,7);
  normas.harm.5 <- rep(0,7);
  normas.harm.6 <- rep(0,7);
  normas.harm.7 <- rep(0,7);
  normas.harm.8 <- rep(0,7);
  normas.harm.9 <- rep(0,7);
  normas.harm.10 <- rep(0,7);
  
  # compute distinct norms of the errors made in the whole graph 
  # use base:: since package Matrix a method called norm
  # 6 is RMSD
  # 7 is MAE (mean absolute error)
  
  normas.harm.1[1] <- base::norm(as.matrix(err.harm.1),type='O');
  normas.harm.1[2] <- base::norm(as.matrix(err.harm.1),type='F');
  normas.harm.1[3] <- base::norm(as.matrix(err.harm.1),type='M');
  normas.harm.1[4] <- cor(est.harm.1,harm.cent,method="spearman");
  normas.harm.1[5] <- cor.test(est.harm.1,harm.cent,method="spearman")$p.value;
  normas.harm.1[6] <- sqrt(mean(err.harm.1^2))
  normas.harm.1[7] <- mean(abs(err.harm.1))
  
  normas.harm.2[1] <- base::norm(as.matrix(err.harm.2),type='O');
  normas.harm.2[2] <- base::norm(as.matrix(err.harm.2),type='F');
  normas.harm.2[3] <- base::norm(as.matrix(err.harm.2),type='M');
  normas.harm.2[4] <- cor(est.harm.2,harm.cent,method="spearman");
  normas.harm.2[5] <- cor.test(est.harm.2,harm.cent,method="spearman")$p.value;
  normas.harm.2[6] <- sqrt(mean(err.harm.2^2))
  normas.harm.2[7] <- mean(abs(err.harm.2))
  
  normas.harm.3[1] <- base::norm(as.matrix(err.harm.3),type='O');
  normas.harm.3[2] <- base::norm(as.matrix(err.harm.3),type='F');
  normas.harm.3[3] <- base::norm(as.matrix(err.harm.3),type='M');
  normas.harm.3[4] <- cor(est.harm.3,harm.cent,method="spearman");
  normas.harm.3[5] <- cor.test(est.harm.3,harm.cent,method="spearman")$p.value;
  normas.harm.3[6] <- sqrt(mean(err.harm.3^2))
  normas.harm.3[7] <- mean(abs(err.harm.3))
  
  normas.harm.4[1] <- base::norm(as.matrix(err.harm.4),type='O');
  normas.harm.4[2] <- base::norm(as.matrix(err.harm.4),type='F');
  normas.harm.4[3] <- base::norm(as.matrix(err.harm.4),type='M');
  normas.harm.4[4] <- cor(est.harm.4,harm.cent,method="spearman");
  normas.harm.4[5] <- cor.test(est.harm.4,harm.cent,method="spearman")$p.value;
  normas.harm.4[6] <- sqrt(mean(err.harm.4^2))
  normas.harm.4[7] <- mean(abs(err.harm.4))
  
  normas.harm.5[1] <- base::norm(as.matrix(err.harm.5),type='O');
  normas.harm.5[2] <- base::norm(as.matrix(err.harm.5),type='F');
  normas.harm.5[3] <- base::norm(as.matrix(err.harm.5),type='M');
  normas.harm.5[4] <- cor(est.harm.5,harm.cent,method="spearman");
  normas.harm.5[5] <- cor.test(est.harm.5,harm.cent,method="spearman")$p.value;
  normas.harm.5[6] <- sqrt(mean(err.harm.5^2))
  normas.harm.5[7] <- mean(abs(err.harm.5))
  
  normas.harm.6[1] <- base::norm(as.matrix(err.harm.6),type='O');
  normas.harm.6[2] <- base::norm(as.matrix(err.harm.6),type='F');
  normas.harm.6[3] <- base::norm(as.matrix(err.harm.6),type='M');
  normas.harm.6[4] <- cor(est.harm.6,harm.cent,method="spearman");
  normas.harm.6[5] <- cor.test(est.harm.6,harm.cent,method="spearman")$p.value;
  normas.harm.6[6] <- sqrt(mean(err.harm.6^2))
  normas.harm.6[7] <- mean(abs(err.harm.6))
  
  normas.harm.7[1] <- base::norm(as.matrix(err.harm.7),type='O');
  normas.harm.7[2] <- base::norm(as.matrix(err.harm.7),type='F');
  normas.harm.7[3] <- base::norm(as.matrix(err.harm.7),type='M');
  normas.harm.7[4] <- cor(est.harm.7,harm.cent,method="spearman");
  normas.harm.7[5] <- cor.test(est.harm.7,harm.cent,method="spearman")$p.value;
  normas.harm.7[6] <- sqrt(mean(err.harm.7^2))
  normas.harm.7[7] <- mean(abs(err.harm.7))
  
  normas.harm.8[1] <- base::norm(as.matrix(err.harm.8),type='O');
  normas.harm.8[2] <- base::norm(as.matrix(err.harm.8),type='F');
  normas.harm.8[3] <- base::norm(as.matrix(err.harm.8),type='M');
  normas.harm.8[4] <- cor(est.harm.8,harm.cent,method="spearman");
  normas.harm.8[5] <- cor.test(est.harm.8,harm.cent,method="spearman")$p.value;
  normas.harm.8[6] <- sqrt(mean(err.harm.8^2))
  normas.harm.8[7] <- mean(abs(err.harm.8))
  
  normas.harm.9[1] <- base::norm(as.matrix(err.harm.9),type='O');
  normas.harm.9[2] <- base::norm(as.matrix(err.harm.9),type='F');
  normas.harm.9[3] <- base::norm(as.matrix(err.harm.9),type='M');
  normas.harm.9[4] <- cor(est.harm.9,harm.cent,method="spearman");
  normas.harm.9[5] <- cor.test(est.harm.9,harm.cent,method="spearman")$p.value;
  normas.harm.9[6] <- sqrt(mean(err.harm.9^2))
  normas.harm.9[7] <- mean(abs(err.harm.9))
  
  normas.harm.10[1] <- base::norm(as.matrix(err.harm.10),type='O');
  normas.harm.10[2] <- base::norm(as.matrix(err.harm.10),type='F');
  normas.harm.10[3] <- base::norm(as.matrix(err.harm.10),type='M');
  normas.harm.10[4] <- cor(est.harm.10,harm.cent,method="spearman");
  normas.harm.10[5] <- cor.test(est.harm.10,harm.cent,method="spearman")$p.value;
  normas.harm.10[6] <- sqrt(mean(err.harm.10^2))
  normas.harm.10[7] <- mean(abs(err.harm.10))
  
  normas.harm <- rbind(normas.harm.1,
                       normas.harm.2,
                       normas.harm.3,
                       normas.harm.4,
                       normas.harm.5,
                       normas.harm.6,
                       normas.harm.7,
                       normas.harm.8,
                       normas.harm.9,
                       normas.harm.10);
 
  data <- list(harm.cent,
               alfa,
               cent_comput_time,
               harm.cuant.1,
               harm.cuant.2,
               harm.cuant.3,
               harm.cuant.4,
               harm.cuant.5,
               harm.cuant.6,
               harm.cuant.7,
               #11
               harm.cuant.8,
               harm.cuant.9,
               harm.cuant.10,
               pistas.harm,
               t_pistas.harm,
               vect_prob,
               cl_instr,
               p.ord.harm.1,
               p.ord.harm.2,
               p.ord.harm.3,
               #21
               p.ord.harm.4,
               p.ord.harm.5,
               p.ord.harm.6,
               p.ord.harm.7,
               p.ord.harm.8,
               p.ord.harm.9,
               p.ord.harm.10,
               est.harm.index.1,
               est.harm.index.2,
               est.harm.index.3,
               #31
               est.harm.index.4,
               est.harm.index.5,
               est.harm.index.6,
               est.harm.index.7,
               est.harm.index.8,
               est.harm.index.9,
               est.harm.index.10,
               est.harm.1,
               est.harm.2,
               est.harm.3,
               #41
               est.harm.4,
               est.harm.5,
               est.harm.6,
               est.harm.7,
               est.harm.8,
               est.harm.9,
               est.harm.10,
               err.harm.1,
               err.harm.2,
               err.harm.3,
               # 51
               err.harm.4,
               err.harm.5,
               err.harm.6,
               err.harm.7,
               err.harm.8,
               err.harm.9,
               err.harm.10,
               normas.harm
  );
  
  return(data);
}

# plot error made by each regression method
# receives proportion of the dataset used as training set, numeric between 0 and 1
# ylimsup_MAE is the upper limit for the y-axis on the MAE plot
graph_err_by_method_4 <- function(file_name,train_size,ylimsup_MAE=NA){
  
  results.simulations <- read.csv(file_name)
  
  # number of networks 
  ncasos <- nrow(results.simulations);
  
  # create data.frame with dataset required to plot
  # group is a factor indicating the method of regression
  df <- data.frame(matrix(ncol = 0, nrow = ncasos * 4))
  
  df[,"group"] <- c(rep("QC",ncasos),
                    rep("L",ncasos),
                    rep("T",ncasos),
                    rep("NN",ncasos));
  
  df[,"L1"] <- c(results.simulations[,2],
                 results.simulations[,3],
                 results.simulations[,6],
                 results.simulations[,8]);
  df[,"L2"] <- c(results.simulations[,2 + 8*1],
                 results.simulations[,3 + 8*1],
                 results.simulations[,6 + 8*1],
                 results.simulations[,8 + 8*1]);
  df[,"Linf"] <- c(results.simulations[,2 + 8*2],
                   results.simulations[,3 + 8*2],
                   results.simulations[,6 + 8*2],
                   results.simulations[,8  + 8*2]);
  df[,"corr"] <- c(results.simulations[,2 + 8*3],
                   results.simulations[,3 + 8*3],
                   results.simulations[,6 + 8*3],
                   results.simulations[,8 + 8*3]);
  df[,"pval"] <- c(results.simulations[,2 + 8*4],
                   results.simulations[,3 + 8*4],
                   results.simulations[,6 + 8*4],
                   results.simulations[,8 + 8*4]);
  df[,"rmsd"] <- c(results.simulations[,2 + 8*5],
                   results.simulations[,3 + 8*5],
                   results.simulations[,6 + 8*5],
                   results.simulations[,8 + 8*5]);
  df[,"mae"] <- c(results.simulations[,2  + 8*6],
                  results.simulations[,3 + 8*6],
                  results.simulations[,6 + 8*6],
                  results.simulations[,8 + 8*6]);
  
  if (is.na(ylimsup_MAE)){
    print(ggplot(df, aes(x=group, y=mae)) + 
            #geom_violin(trim=TRUE) + 
            geom_boxplot(width=.1, aes(fill=group),outlier.size=1.5, outlier.shape=21) +
            stat_summary(fun=mean, geom="point", fill="white", shape=23, size=2.5)+
            labs(x="Regression method", y=paste("Distribution across ",toString(ncasos)," graphs", sep = "")#,
                 #title=paste("Estimation error - MAE - Training Size = ",toString(train_size * 100)," %",sep="")
            )+
            theme(plot.title = element_text(size=15, face="bold",margin = margin(10, 0, 10, 0),
                                            hjust=0.5))
    );
  }
  else{
    print(ggplot(df, aes(x=group, y=mae)) + 
            #geom_violin(trim=TRUE) + 
            geom_boxplot(width=.1, aes(fill=group),outlier.size=1.5, outlier.shape=21) +
            coord_cartesian(ylim = c(0,ylimsup_MAE)) +
            stat_summary(fun=mean, geom="point", fill="white", shape=23, size=2.5)+
            labs(x="Regression method", y=paste("Distribution across ",toString(ncasos)," graphs", sep = "")#,
                 #title=paste("Estimation error - MAE - Training Size = ",toString(train_size * 100)," %",sep="")
            )+
            theme(plot.title = element_text(size=15, face="bold",margin = margin(10, 0, 10, 0),
                                            hjust=0.5))
    );
  }
  
  return(df);
}

# plot error in function of the dataset proportion used as training set
# receives name of .csv file storing results 
graph_err_by_siz_3 <- function(file_name){
  
  results.simulations <- read.csv(file_name)
  
  # number of networks
  ncasos <- nrow(results.simulations);
  
  # build data.frame with required data to plot
  # group is a factor indicating the proportion of the dataset used as training set
  df <- data.frame(matrix(ncol = 0, nrow = ncasos * 10))
  
  df[,"group"] <- c(rep("1.0",ncasos),
                    rep("0.9",ncasos),
                    rep("0.8",ncasos),
                    rep("0.7",ncasos),
                    rep("0.6",ncasos),
                    rep("0.5",ncasos),
                    rep("0.4",ncasos),
                    rep("0.3",ncasos),
                    rep("0.2",ncasos),
                    rep("0.1",ncasos));
  
  df[,"L1"] <- c(results.simulations[,2],
                 results.simulations[,3],
                 results.simulations[,4],
                 results.simulations[,5],
                 results.simulations[,6],
                 results.simulations[,7],
                 results.simulations[,8],
                 results.simulations[,9],
                 results.simulations[,10],
                 results.simulations[,11]);
  
  df[,"L2"] <- c(results.simulations[,2+ 1 * 10],
                 results.simulations[,3+ 1 * 10],
                 results.simulations[,4+ 1 * 10],
                 results.simulations[,5+ 1 * 10],
                 results.simulations[,6+ 1 * 10],
                 results.simulations[,7+ 1 * 10],
                 results.simulations[,8+ 1 * 10],
                 results.simulations[,9+ 1 * 10],
                 results.simulations[,10+ 1 * 10],
                 results.simulations[,11+ 1 * 10]);
  
  df[,"Linf"] <- c(results.simulations[,2 + 2 * 10],
                   results.simulations[,3 + 2 * 10],
                   results.simulations[,4 + 2 * 10],
                   results.simulations[,5 + 2 * 10],
                   results.simulations[,6 + 2 * 10],
                   results.simulations[,7 + 2 * 10],
                   results.simulations[,8 + 2 * 10],
                   results.simulations[,9 + 2 * 10],
                   results.simulations[,10 + 2 * 10],
                   results.simulations[,11 + 2 * 10]);
  
  df[,"corr"] <- c(results.simulations[,2 + 3 * 10],
                   results.simulations[,3 + 3 * 10],
                   results.simulations[,4 + 3 * 10],
                   results.simulations[,5 + 3 * 10],
                   results.simulations[,6 + 3 * 10],
                   results.simulations[,7 + 3 * 10],
                   results.simulations[,8 + 3 * 10],
                   results.simulations[,9 + 3 * 10],
                   results.simulations[,10 + 3 * 10],
                   results.simulations[,11 + 3 * 10]);
  
  df[,"pval"] <- c(results.simulations[,2 + 4 * 10],
                   results.simulations[,3 + 4 * 10],
                   results.simulations[,4 + 4 * 10],
                   results.simulations[,5 + 4 * 10],
                   results.simulations[,6 + 4 * 10],
                   results.simulations[,7 + 4 * 10],
                   results.simulations[,8 + 4 * 10],
                   results.simulations[,9 + 4 * 10],
                   results.simulations[,10 + 4 * 10],
                   results.simulations[,11 + 4 * 10]);
  
  df[,"rmsd"] <- c(results.simulations[,2 + 5 * 10],
                   results.simulations[,3 + 5 * 10],
                   results.simulations[,4 + 5 * 10],
                   results.simulations[,5 + 5 * 10],
                   results.simulations[,6 + 5 * 10],
                   results.simulations[,7 + 5 * 10],
                   results.simulations[,8 + 5 * 10],
                   results.simulations[,9 + 5 * 10],
                   results.simulations[,10 + 5 * 10],
                   results.simulations[,11 + 5 * 10]);
  
  df[,"mae"] <- c(results.simulations[,2 + 6 * 10],
                  results.simulations[,3 + 6 * 10],
                  results.simulations[,4 + 6 * 10],
                  results.simulations[,5 + 6 * 10],
                  results.simulations[,6 + 6 * 10],
                  results.simulations[,7 + 6 * 10],
                  results.simulations[,8 + 6 * 10],
                  results.simulations[,9 + 6 * 10],
                  results.simulations[,10 + 6 * 10],
                  results.simulations[,11 + 6 * 10]);
  
  print(ggplot(df, aes(x=group, y=mae)) + 
          #geom_violin(trim=TRUE) + 
          geom_boxplot(width=.1, aes(fill=group),outlier.size=1.5, outlier.shape=21) +
          stat_summary(fun=mean, geom="point", fill="white", shape=23, size=2.5)+
          labs(x="Proportion of nodes employed to estimate alpha", 
               y=paste("Distribution across ",toString(ncasos)," graphs", sep = "")#,
               #title="Estimation error - MAE"
          )+
          theme(plot.title = element_text(size=15, face="bold",margin = margin(10, 0, 10, 0),
                                          hjust=0.5)))
  
  return(df);
}

# simulations to check that harmonic centrality in BA graphs holds the required assumptions
# n_iter is the number of networks
assumption_validation <- function(n_iter,n_it_boot){
  
  # ncol is the number of statistics to be recorded
  df <- data.frame(matrix(ncol = 8, nrow = 0))
  
  n_cores <- parallel::detectCores();
  
  for(i in 1:n_iter){
    
    gr <- sample_pa(10000,power=1);
    har <- harmonic(gr);
    deg <- degree(gr);
    
    # nodes with 0 centrality are removed
    har_0 <- har[which(har!=0)];
    # create object poweRlaw
    m_m_0 <- conpl$new(har_0)
    
    # estimate x_min and alpha by minimizing KS
    est <- estimate_xmin(m_m_0)
    # store x_min and alpha
    df[i,1] =est$xmin;
    df[i,2] =est$pars;
    df[i,3] = est$gof;
    # bootstrap to get significance of the gof 
    bs_p = bootstrap_p(m_m_0, no_of_sims=n_it_boot, threads=n_cores)
    df[i,4] =bs_p$p;
    
    # set x_min to 1
    m_m_0$setXmin(1);
    # recalculate alpha
    est2 <- estimate_pars(m_m_0)
    df[i,5] =est2$pars;
    
    # compute correlation and significance
    index_no_0 <- union(which(deg!=0),which(har!=0));
    har_1 <- log(har[index_no_0]);
    deg_1 <- log(deg[index_no_0]);
    df[i,6] = cor(har_1,deg_1,method="spearman");
    df[i,7] = cor.test(har_1,deg_1,method="spearman")$p.value;
    
    df[i,8] = min(har_0)
  }
  
  return(df)
}

# obtains the probability quantiles associated to the values in vector central
# producing equi-spaced variability of central in a number of n_cuant values
obtener_cuantiles <- function(central,n_cuant){
  mi <- min(central)
  ma <- max(central)
  max_var <- seq(mi,ma,by=(ma-mi)/(n_cuant+1))
  
  sorted_central <- sort(central)
  lsc <- length(sorted_central)
  
  # values in central nearest to the values in max_var
  max_var_cent_real <- sapply(max_var,function(x) 
  {binary_search(x,sorted_central)})
  
  # find position of these values in the sorted version of central, to get the quantiles
  cent.r <- match(max_var_cent_real,sorted_central)
  
  return(list((cent.r/lsc)[1:n_cuant + 1],max_var_cent_real))
}

# looks for the nearest value in vector to val
# by means of binary search, vector is sorted in ascending order
binary_search <- function(val, vector){
  l <- length(vector);
  
  # base case: vector of length 2
  if(l==2){
    if(abs(vector[1]-val)<abs(vector[2]-val))
      return(vector[1])
    else{return(vector[2])}
  }
  
  l2 <- round(l/2)
  
  if(vector[l2]<val)
    return(binary_search(val,vector[l2:l]))
  else return(binary_search(val,vector[1:l2]))
}

# compute distinct error norms
compute_error_norms <- function(estim,real){
  
  err <- estim - real;
  
  normas <- rep(0,7);
  
  # several norms of the error made in the whole graph
  # use base:: since library Matrix has another method named norm
  # 1 L1
  # 2 L2
  # 3 Linf
  # 6 es RMSD
  # 7 es MAE (mean absolute error)
  
  normas[1] <- base::norm(as.matrix(err),type='O');
  normas[2] <- base::norm(as.matrix(err),type='F');
  normas[3] <- base::norm(as.matrix(err),type='M');
  normas[4] <- cor(estim,real,method="spearman");
  normas[5] <- cor.test(estim,real,method="spearman")$p.value;
  normas[6] <- sqrt(mean(err^2))
  normas[7] <- mean(abs(err))
  
  return(normas);
  
}

# results required to produce Fig 1
MAE_dist_train_sizes <- function(){
  harm.cent <- harmonic(sample_pa(10000,power=1))
  harm.cuant <- obtener_cuantiles(log(harm.cent+1),8);
  pv8 <- harm.cuant[[1]]
  
  resultados.simulaciones <- get_error_by_prop(1000,pv8)
  write.csv(resultados.simulaciones,file="resultados_simulaciones8.csv");
}

# code to produce Fig 1
plot_MAE_dist_train_sizes <- function(){
  res_prop_8 <- graph_err_by_siz_3("resultados_simulaciones8.csv")
}

# results required to produce Fig 2, 3
produce_comparison_methods <- function() {
  harm.cent <- harmonic(sample_pa(10000,power=1))
  harm.cuant <- obtener_cuantiles(log(harm.cent+1),8);
  pv8 <- harm.cuant[[1]]
  
  resultados.simulaciones.metodos.1 <- get_error_by_method_3(1000,pv8,1);
  write.csv(resultados.simulaciones.metodos.1,file="resultados_metodos8_1.csv");
  
  resultados.simulaciones.metodos.0.1 <- get_error_by_method_3(1000,pv8,0.1);
  write.csv(resultados.simulaciones.metodos.0.1,file="resultados_metodos8_0.1.csv");
}

# code to produce Fig 2, 3
plot_comparison_methods <- function(){
  graph_err_by_method_4("resultados_metodos8_1.csv",1,8);
  graph_err_by_method_4("resultados_metodos8_0.1.csv",0.1,20);
}

get_statistics_Table_I_II <- function(){
  df <- read.csv("valid_assump.csv")
  df[,10] <- df[,3] - df[,6];
  apply(df,2,summary);
}

# filename = "resultados_metodos8_0.1.csv" or 
# filename = "resultados_metodos8_1.csv"
get_statistics_Table_III_IV <- function(filename) {
  results.simulations <- read.csv(filename)
  elaps_time <- data.frame(matrix(nrow=1000,ncol=0));
  elaps_time[,"L"] <- results.simulations[,3 +0*5 + 8*7 +1]
  # times include construction of clues, medians computation, and QuickCent time
  elaps_time[,"QC"] <- results.simulations[,3 +6*5+ 8*7 +1] + results.simulations[,3 +5*5+ 8*7 +1] + results.simulations[,3 +4*5+ 8*7 +1]
  elaps_time[,"T"] <- results.simulations[,3 +5*1+ 8*7 +1]
  elaps_time[,"NN"] <- results.simulations[,3 +2*5+ 8*7+1]
  return(elaps_time)
}
# and then
#format(apply(elaps_time,2,mean)*1000,digits=2);
#format(apply(elaps_time,2,sd)*1000,digits=2);