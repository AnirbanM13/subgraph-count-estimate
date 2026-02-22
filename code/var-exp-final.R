####################### Theta ##################################################

theta <- function(PI, alpha){
  lambda <- tabulate(alpha)/length(alpha)
  K <- max(alpha)
  
  theta.1 <- 0
  theta.2 <- 0
  
  for (u in 1:K) {
    for (v in 1:K) {
      for (w in 1:K) {
        theta.1 <- PI[u,v]*PI[v,w]*lambda[u]*lambda[v]*lambda[w] + theta.1
        theta.2 <- PI[u,v]*PI[v,w]*PI[w,u]*lambda[u]*lambda[v]*lambda[w] + theta.2
      }
    }
  }
  
  return(c(theta.1, theta.2))
  
}

##################### Induced Network formation ###############################
###############################################################################



# Edge 
var_Edge_ind <- function(N, alpha, Pi, p){
  K <- max(alpha)
  lambda <- tabulate(alpha)/length(alpha)
  
  c1 <- (p^3)*(1-p)
  c2 <- (p^2)*(1-p^2)
  
  sum1 <- 0
  sum2 <- 0
  for (u in 1:K) {
    for (v in 1:K) {
      sum2 <- Pi[u,v]*lambda[u]*lambda[v] + sum2
      for (w in 1:K) {
        sum1 <- Pi[u,v]*Pi[u,w]*lambda[u]*lambda[v]*lambda[w] + sum1
      }
    }
  }
  
  total <- 4*N^3*c1*sum1 + N^2*c2*sum2
  
  return(total)
}


# Wedge
var_wedge_ind <- function(N,alpha,Pi,p){
  
  K <- max(alpha)
  la <- tabulate(alpha)/length(alpha)
  
  sum1 <- 0
  sum2 <- 0
  sum3 <- 0
  
  for (u in 1:K) {
    for (v in 1:K) {
      for (w in 1:K) {
        sum1 <- p^3 * (1 - p^3) * Pi[u,v]*Pi[u,w] * (2  + 4 * Pi[v,w]) * la[u]  * la[v] * la[w] + sum1
        for (x in 1:K) {
          sum2 <- p^4 * (1 - p^2) *Pi[u,v] * Pi[u,w] * (3 *  Pi[u,x]  
                                                        + 3 * Pi[v,x] +
                                                          24  * Pi[v,x] * Pi[u,x] +
                                                          6 * Pi[v,x] * Pi[w,x]) * la[u] * la[v] * la[w] * la[x] + sum2
          for (y in 1:K) {
            sum3 <-p^5 * (1-p) * Pi[u,v] * Pi[u,w] * (Pi[u,x] * Pi[u,y] + 
                                                        4  * Pi[u,x] * Pi[x,y] +
                                                        4  * Pi[v,x] * Pi[x,y]) * la[u] * la[v] * la[w] * la[x] * la[y] + sum3
          }
        }
      }
    }
  }
  
  
  t <- N^3 * sum1 + N^4 * sum2 + N^5 * sum3
  return(t)
}


# Triangle
var_Triangle_ind <- function(N, alpha,Pi,p){
  K <- max(alpha)
  lambda <- tabulate(alpha)/length(alpha)
  
  c3 <- p^3*(1-p^3)
  c2 <- p^4*(1-p^2)
  c1 <- p^5*(1-p)
  
  sum1 <- 0
  sum2 <- 0
  sum3 <- 0
  
  for (u in 1:K) {
    for (v in 1:K) {
      for (w in 1:K) {
        sum3 <- Pi[u,v]*Pi[v,w]*Pi[w,u]*lambda[u]*lambda[v]*lambda[w] + sum3
        for (x in 1:K) {
          sum2 <- Pi[u,v]*Pi[v,w]*Pi[w,u]*Pi[v,x]*Pi[x,u]*lambda[u]*lambda[v]*lambda[w]*lambda[x] + sum2
          for (y in 1:K) {
            sum1 <- (Pi[u,v]*Pi[v,w]*Pi[w,u])*(Pi[u,x]*Pi[x,y]*Pi[y,u])*lambda[u]*lambda[v]*lambda[w]*lambda[x]*lambda[y] + sum1
          }
        }
      }
    }
  }
  
  total <- 9*N^5*c1*sum1 + 18*N^4*c2*sum2 + 6*N^3*c3*sum3 
  return(total)
}

# Covariance between wedge and triangle
cov_wedge_tri_ind <- function(N,alpha,Pi,p){
  
  K <- max(alpha)
  la <- tabulate(alpha)/length(alpha)
  
  sum1 <- 0
  sum2 <- 0
  sum3 <- 0
  
  for (u in 1:K) {
    for (v in 1:K) {
      for (w in 1:K) {
        sum1 <- 6 * p^3 * (1 - p^3) * Pi[u,v] * Pi[v,w] * Pi[w,u] * la[u] * la[v] * la[w] + sum1
        for (x in 1:K) {
          sum2 <- p^4 * (1 - p^2) * Pi[u,v] * Pi[u,w] * (12  * Pi[u,x] * Pi[v,x] +
                                                           6 * Pi[v,w] * Pi[v,x] * Pi[w,x]) * la[u] * la[v] * la[w] * la[x] + sum2
          for (y in 1:K) {
            sum3 <- p^5 * (1 - p) * Pi[u,v] * Pi[u,w] * (3  * Pi[u,x] * Pi[x,y] * Pi[y,u] +
                                                           6 *  Pi[v,x] * Pi[x,y] * Pi[y,v]) * la[u] * la[v] * la[w] *la[x] * la[y] + sum3 
          }
        }
      }
    }
  }
  
  t <- N^3 * sum1 + N^4 * sum2 + N^5 * sum3
  return(t)
}


# Clustering Coefficient
var_cc_ind <- function(N,alpha,Pi,p, beta){
  
  
  v_w <- var_wedge_ind(N,alpha,Pi,p)
  v_t <- var_Triangle_ind(N,alpha,Pi,p)
  cov_wt <- cov_wedge_tri_ind(N,alpha,Pi,p)
  
  theta1 <- N^(2*beta) * theta(Pi, alpha)[1]
  theta2 <- N^(3*beta) * theta(Pi, alpha)[2]
  
  
  t11 <- ((theta2)^2)/((theta1)^4 * N^(6 - 4*beta) * p^6)
  t22 <- 1/((theta1)^2 * N^(6 - 6 * beta) * p^6)
  t12 <- ((theta2))/((theta1)^3 * N^(6 - 5*beta) * p^6)
  
  v_cc <- t11 * v_w + t22 * v_t - 2 * t12 * cov_wt
  
  return(v_cc)
}


################################################################################
########## Ego-centric Network formation #######################################
################################################################################



# Edge 
var_Edge_ego <- function(N,alpha, Pi, p){
  K <- max(alpha)
  lambda <- tabulate(alpha)/length(alpha)
  
  c2 <- p*(2-5*p+4*p^2-p^3)
  c1 <- p*(1-p)^3
  
  sum1 <- 0
  sum2 <- 0
  for (u in 1:K) {
    for (v in 1:K) {
      sum2 <- Pi[u,v]*lambda[u]*lambda[v] + sum2
      for (w in 1:K) {
        sum1 <- Pi[u,v]*Pi[u,w]*lambda[u]*lambda[v]*lambda[w] + sum1
      }
    }
  }
  
  total <- 4*N^3*c1*sum1 + N^2*c2*sum2
  
  return(total)
}



# Wedge
var_wedge_ego <- function(N,alpha,Pi,p){
  
  K <- max(alpha)
  la <- tabulate(alpha)/length(alpha)
  
  sum1 <- 0
  sum2 <- 0
  sum3 <- 0
  
  for (u in 1:K) {
    for (v in 1:K) {
      for (w in 1:K) {
        sum1 <- Pi[u,v]*Pi[u,w] * (2 * (p + (1-p)*p^2)*(1 - (p + (1-p)*p^2)) 
                                   + 4 * p^2 * (1-p)^2 * (2-p^2) * Pi[v,w])  * la[u] * la[v] * la[w] + sum1
        for (x in 1:K) {
          sum2 <- Pi[u,v] * Pi[u,w] * (3 * p * (1-p)^2 * (1 + p - p^3) * Pi[u,x]  
                                       + 3 * p^2 * (1 - p) * (2 - p^2) * Pi[v,x] +
                                         24 * p^2 * (1 - p)^3 * Pi[v,x] * Pi[u,x] +
                                         6 * p^2 * (1 - p)^3 * (1 + p) * Pi[v,x] * Pi[w,x]) * la[u] * la[v] * la[w] * la[x] + sum2
          for (y in 1:K) {
            sum3 <- Pi[u,v] * Pi[u,w] * (p * (1 - p)^3* (1+p)^2 * Pi[u,x] * Pi[u,y] + 
                                           4 * p^2 * (1 - p)^3 * (1 + p) * Pi[u,x] * Pi[x,y] +
                                           4 * p^3 * (1 - p)^3 * Pi[v,x] * Pi[x,y]) * la[u] * la[v] * la[w] * la[x] * la[y] + sum3
          }
        }
      }
    }
  }
  
  
  t <- N^3 * sum1 + N^4 * sum2 + N^5 * sum3
  return(t)
}



# Triangle
var_Triangle_ego <- function(N,alpha,Pi,p){
  K <- max(alpha)
  lambda <- tabulate(alpha)/length(alpha)
  
  c3 <- (3*p^2 - 2*p^3)*(1-(3*p^2 - 2*p^3))
  c2 <- p^2 * (1-p) * (1 - p -8*p^2 +4*p^3)
  c1 <- 4*p^3*(1-p)^3
  
  sum1 <- 0
  sum2 <- 0
  sum3 <- 0
  
  for (u in 1:K) {
    for (v in 1:K) {
      for (w in 1:K) {
        sum3 <- Pi[u,v]*Pi[v,w]*Pi[w,u]*lambda[u]*lambda[v]*lambda[w] + sum3
        for (x in 1:K) {
          sum2 <- Pi[u,v]*Pi[v,w]*Pi[w,u]*Pi[v,x]*Pi[x,u]*lambda[u]*lambda[v]*lambda[w]*lambda[x] + sum2
          for (y in 1:K) {
            sum1 <- Pi[u,v]*Pi[v,w]*Pi[w,u]*Pi[u,x]*Pi[x,y]*Pi[y,u]*lambda[u]*lambda[v]*lambda[w]*lambda[x]*lambda[y] + sum1
          }
        }
      }
    }
  }
  
  total <- 9*N^5*c1*sum1 + 18*N^4*c2*sum2 + 6*N^3*c3*sum3 
  return(total)
}



# Covariance between wedge and triangle

cov_wedge_tri_ego <- function(N,alpha,Pi,p){
  
  K <- max(alpha)
  la <- tabulate(alpha)/length(alpha)
  
  sum1 <- 0
  sum2 <- 0
  sum3 <- 0
  
  for (u in 1:K) {
    for (v in 1:K) {
      for (w in 1:K) {
        sum1 <- 6 * p^2 *(1 - p) * (3- 3*p -3*p^2 +2*p^3) * Pi[u,v] * Pi[v,w] * Pi[w,u] * la[u] * la[v] * la[w] + sum1
        for (x in 1:K) {
          sum2 <- Pi[u,v] * Pi[u,w] * (12 * p^2 * (1 - p) * (2 -p -3*p^2 +2*p^3) * Pi[u,x] * Pi[v,x] +
                                         6 * p^2 * (1 - p)^2 * (1+ p-2*p^2) *Pi[v,w]* Pi[v,x] * Pi[w,x]) * la[u] * la[v] * la[w] * la[x] + sum2
          for (y in 1:K) {
            sum3 <- Pi[u,v] * Pi[u,w] * ((3 *2 * p^2 * (1 - p)^3 * (1+ p) * Pi[u,x] * Pi[x,y] * Pi[y,u]) +
                                           (6 * 2* p^3 * (1 - p)^3 * Pi[v,x] * Pi[x,y] * Pi[y,v])) * la[u] * la[v] * la[w] *la[x] * la[y] + sum3 
          }
        }
      }
    }
  }
  
  t <- N^3 * sum1 + N^4 * sum2 + N^5 * sum3
  return(t)
}


# Clustering coefficient

var_cc_ego <- function(N,alpha,Pi,p, beta){
  
  v_w <- var_wedge_ego(N,alpha,Pi,p)
  v_t <- var_Triangle_ego(N,alpha,Pi,p)
  cov_wt <- cov_wedge_tri_ego(N,alpha,Pi,p)
  
  theta1 <- N^(2*beta) * theta(Pi, alpha)[1]
  theta2 <- N^(3*beta) * theta(Pi, alpha)[2]
  
  
  t11 <- ((theta2)^2)/((theta1)^4 * N^(6 - 4*beta) * (p + (1-p) * p^2)^2)
  t22 <- 1/((theta1)^2 * N^(6 - 6 * beta) * (3*p^2  - 2*p^3)^2)
  t12 <- ((theta2))/((theta1)^3 * N^(6 - 5*beta) * (p + (1-p) * p^2) * (3*p^2  - 2*p^3))
  
  v_cc <- t11 * v_w + t22 * v_t - 2 * t12 * cov_wt
  
  return(v_cc)
}













