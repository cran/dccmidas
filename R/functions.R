#' Power of a matrix
#'
#' Reports the power of a matrix.
#' @param x A square matrix
#' @param p The power of interest
#' @return The resulting matrix x^(p)
#' @importFrom Rdpack reprompt
#' @keywords internal
#' @export

"%^%" <- function(x, p) 
   with(eigen(x), vectors %*% (values^p* t(vectors)))

#' DCC-MIDAS log-likelihood (second step)
#'
#' Obtains the log-likelihood of the DCC models in the second step.
#' For details, see \insertCite{colacito2011component;textual}{dccmidas} and \insertCite{dcc_engle_2002;textual}{dccmidas}.
#' @param param Vector of starting values. 
#' @param res Array of standardized daily returns, coming from the first step estimation.
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively.
#' @param N_c Number of (lagged) realizations to use for the standarized residuals forming the long-run correlation.
#' @param K_c Number of (lagged) realizations to use for the long-run correlation.
#' @return The resulting vector is the log-likelihood value for each \eqn{t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @import rumidas
#' @references
#' \insertAllCited{} 
#' @keywords internal
#' @export

dccmidas_loglik<-function(param,res,lag_fun="Beta",N_c,K_c){

a<-param[1]
b<-param[2]
w1<- ifelse(lag_fun=="Beta",1,0)
w2<-param[3]

##### index

Num_col<-dim(res)[1]		# number of assets
TT<-dim(res)[3]			# number of daily observations

##### matrices and vectors

C_t<-array(0,dim=c(Num_col,Num_col,TT))
V_t<-array(0,dim=c(Num_col,Num_col,TT))
Prod_eps_t<-array(0,dim=c(Num_col,Num_col,TT))
Q_t<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
R_t<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
log_det_R_t<-rep(0,TT)
R_t_solved<-array(1,dim=c(Num_col,Num_col,TT))
Eps_t_cross_prod<-rep(0,TT)
Eps_t_R_t_Eps_t<-rep(0,TT)

##### first cycle

for(tt in (N_c+1):TT){
V_t[,,tt]<-rowSums(res[,1,tt:(tt-N_c)]*(res[,1,tt:(tt-N_c)]))*diag(Num_col)
Prod_eps_t[,,tt]<-res[,,tt:(tt-N_c)]%*%t(res[,,tt:(tt-N_c)])
V_t_0.5<-Inv(sqrt(V_t[,,tt]))
C_t[,,tt]<-V_t_0.5%*%Prod_eps_t[,,tt]%*%V_t_0.5
}

#### R_t_bar

weight_fun<-ifelse(lag_fun=="Beta",rumidas::beta_function,rumidas::exp_almon)

betas<-c(rev(weight_fun(1:(K_c+1),(K_c+1),w1,w2))[2:(K_c+1)],0)
R_t_bar<-array(1,dim=c(Num_col,Num_col,TT))

matrix_id<-matrix(1:Num_col^2,ncol=Num_col)
matrix_id_2<-which(matrix_id==1:Num_col^2, arr.ind=TRUE)

for(i in 1:nrow(matrix_id_2)){
R_t_bar[matrix_id_2[i,1],matrix_id_2[i,2],]                      <- suppressWarnings(
roll::roll_sum(C_t[matrix_id_2[i,1],matrix_id_2[i,2],], c(K_c+1),weights = betas)
) 
}

R_t_bar[,,1:K_c]<-diag(Num_col)

################# likelihood

ll<-rep(0,TT)

for(tt in (K_c+1):TT){
Q_t[,,tt]<-(1-a-b)*R_t_bar[,,tt] + a*res[,,tt-1]%*%t(res[,,tt-1]) + b*Q_t[,,tt-1]
Q_t_star<-Inv(sqrt(diag(diag(Q_t[,,tt]))))
R_t[,,tt]<-Q_t_star%*%Q_t[,,tt]%*%Q_t_star
log_det_R_t[tt]<-log(Det(R_t[,,tt]))
R_t_solved[,,tt]<-Inv(R_t[,,tt])
#Eps_t_cross_prod[tt]<-rbind(res[,,tt])%*%cbind(res[,,tt])
Eps_t_R_t_Eps_t[tt]<-rbind(res[,,tt])%*%R_t_solved[,,tt]%*%cbind(res[,,tt])
}

ll<- - (log_det_R_t+Eps_t_R_t_Eps_t)

#sum(ll)

return(ll)

}

#' Obtains the matrix H_t, R_t and long-run correlations, under the DCC-MIDAS model
#'
#' Obtains the matrix H_t, R_t and long-run correlations, under the DCC-MIDAS model
#' For details, see \insertCite{colacito2011component;textual}{dccmidas} and \insertCite{dcc_engle_2002;textual}{dccmidas}.
#' @param est_param Vector of estimated values
#' @param res Array of standardized daily returns, coming from the first step estimation
#' @param Dt Matrix of conditional standard deviations (coming from the first step)
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively
#' @param N_c Number of (lagged) realizations to use for the standarized residuals forming the long-run correlation
#' @param K_c Number of (lagged) realizations to use for the long-run correlation
#' @return A list with the \eqn{H_t}, \eqn{R_t} and long-run correlaton matrices, for each \eqn{t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @import rumidas
#' @references
#' \insertAllCited{} 
#' @keywords internal
#' @export


dccmidas_mat_est<-function(est_param,res,Dt,lag_fun="Beta",N_c,K_c){

a<-est_param[1]
b<-est_param[2]
w1<-ifelse(lag_fun=="Beta",1,0)
w2<-est_param[3]

##### index

Num_col<-dim(res)[1]		# number of assets
TT<-dim(res)[3]			# number of daily observations

##### matrices and vectors

C_t<-array(0,dim=c(Num_col,Num_col,TT))
V_t<-array(0,dim=c(Num_col,Num_col,TT))
Prod_eps_t<-array(0,dim=c(Num_col,Num_col,TT))
Q_t<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
R_t<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
log_det_R_t<-rep(0,TT)
R_t_solved<-array(1,dim=c(Num_col,Num_col,TT))
Eps_t_cross_prod<-rep(0,TT)
Eps_t_R_t_Eps_t<-rep(0,TT)

S<-stats::cov(t(apply(res, 3L, c)))
H_t<-array(S,dim=c(Num_col,Num_col,TT))

##### first cycle

for(tt in (N_c+1):TT){
V_t[,,tt]<-rowSums(res[,1,tt:(tt-N_c)]*(res[,1,tt:(tt-N_c)]))*diag(Num_col)
Prod_eps_t[,,tt]<-res[,,tt:(tt-N_c)]%*%t(res[,,tt:(tt-N_c)])
V_t_0.5<-Inv(sqrt(V_t[,,tt]))
C_t[,,tt]<-V_t_0.5%*%Prod_eps_t[,,tt]%*%V_t_0.5
}

#### R_t_bar

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas<-c(rev(weight_fun(1:(K_c+1),(K_c+1),w1,w2))[2:(K_c+1)],0)
R_t_bar<-array(1,dim=c(Num_col,Num_col,TT))

matrix_id<-matrix(1:Num_col^2,ncol=Num_col)
matrix_id_2<-which(matrix_id==1:Num_col^2, arr.ind=TRUE)

for(i in 1:nrow(matrix_id_2)){
R_t_bar[matrix_id_2[i,1],matrix_id_2[i,2],]                      <- suppressWarnings(
roll_sum(C_t[matrix_id_2[i,1],matrix_id_2[i,2],], c(K_c+1),weights = betas)
) 
}

R_t_bar[,,1:K_c]<-diag(rep(1,Num_col))

################# H_t, R_t and R_t_bar

for(tt in (K_c+1):TT){
Q_t[,,tt]<-(1-a-b)*R_t_bar[,,tt] + a*res[,,tt-1]%*%t(res[,,tt-1]) + b*Q_t[,,tt-1]
Q_t_star<-Inv(sqrt(diag(diag(Q_t[,,tt]))))
R_t[,,tt]<-Q_t_star%*%Q_t[,,tt]%*%Q_t_star
H_t[,,tt]<-Dt[,,tt]%*%R_t[,,tt]%*%Dt[,,tt]
}

#H_t_m<-apply(H_t[,,(K_c+1):TT],c(1,2),mean)
#H_t[,,1:K_c]<-H_t_m

results<-list(
"H_t"=H_t,
"R_t"=R_t,
"R_t_bar"=R_t_bar
)

return(results)

}

#' RiskMetrics model
#'
#' Obtains the matrix H_t, under the RiskMetrics model.
#' @param r_t List of daily returns
#' @param lambda **optional** Decay parameter. Default to 0.94
#' @return A list with the \eqn{H_t} matrix, for each \eqn{t}.
#' @importFrom Rdpack reprompt
#' \donttest{
#' require(xts)
#' # close to close daily log-returns
#' r_t_s<-diff(log(sp500['2010/2019'][,3]))
#' r_t_s[1]<-0
#' r_t_n<-diff(log(nasdaq['2010/2019'][,3]))
#' r_t_n[1]<-0
#' r_t_f<-diff(log(ftse100['2010/2019'][,3]))
#' r_t_f[1]<-0
#' db_m<-merge.xts(r_t_s,r_t_n,r_t_f)
#' db_m<-db_m[complete.cases(db_m),]
#' colnames(db_m)<-c("S&P500","NASDAQ","FTSE100")
#' # list of returns
#' r_t<-list(db_m[,1],db_m[,2],db_m[,3])
#' RM<-riskmetrics_mat(r_t)
#' }
#' @export


riskmetrics_mat<-function(r_t,lambda=0.94){

################################### checks

cond_r_t<- class(r_t)

if(cond_r_t != "list") { stop(
cat("#Warning:\n Parameter 'r_t' must be a list object. Please provide it in the correct form \n")
)}

if(lambda<0|lambda>1) { stop(
cat("#Warning:\n Parameter 'lambda' must be in the interval 0-1 \n")
)}


################################### merge (eventually)

k<-length(r_t)

len<-rep(NA,k)

for(i in 1:k){
len[i]<-length(r_t[[i]])
}

db<- do.call(xts::merge.xts,r_t)
db<-db[stats::complete.cases(db),] 

for(i in 1:k){
colnames(db)[i]<-colnames(r_t[[i]])
}

TT<-nrow(db)

db<-zoo::coredata(db)

##### matrices and vectors

S<-stats::cov(db)
H_t<-array(S,dim=c(k,k,TT))
Prod_r_t<-array(0,dim=c(k,k,TT))

##### cycle

for(tt in (2):TT){
Prod_r_t[,,tt-1]<-db[(tt-1),]%*%t(db[(tt-1),])
H_t[,,tt]<-lambda*Prod_r_t[,,tt-1] + (1-lambda)*H_t[,,tt-1]
}

results<-list(
"H_t"=H_t
)

return(results)

}


#' Moving Covariance model
#'
#' Obtains the matrix H_t, under the Moving Covariance model.
#' @param r_t List of daily returns
#' @param V  Length of the rolling window adopted. By default, V is 22
#' @return A list with the \eqn{H_t} matrix, for each \eqn{t}.
#' @importFrom Rdpack reprompt
#' @examples
#' \donttest{
#' require(xts)
#' # close to close daily log-returns
#' r_t_s<-diff(log(sp500['2010/2019'][,3]))
#' r_t_s[1]<-0
#' r_t_n<-diff(log(nasdaq['2010/2019'][,3]))
#' r_t_n[1]<-0
#' r_t_f<-diff(log(ftse100['2010/2019'][,3]))
#' r_t_f[1]<-0
#' db_m<-merge.xts(r_t_s,r_t_n,r_t_f)
#' db_m<-db_m[complete.cases(db_m),]
#' colnames(db_m)<-c("S&P500","NASDAQ","FTSE100")
#' # list of returns
#' r_t<-list(db_m[,1],db_m[,2],db_m[,3])
#' MC<-moving_cov(r_t,V=60)
#' }
#' @export

moving_cov<-function(r_t,V=22){

################################### checks

cond_r_t<- class(r_t)

if(cond_r_t != "list") { stop(
cat("#Warning:\n Parameter 'r_t' must be a list object. Please provide it in the correct form \n")
)}

if(V<length(r_t))
 { stop(
cat("#Warning:\n Parameter 'V' must be greater than the number of assets in r_t \n")
)}


################################### merge (eventually)

k<-length(r_t)

len<-rep(NA,k)

for(i in 1:k){
len[i]<-length(r_t[[i]])
}


db<-do.call(xts::merge.xts,r_t)
db<-db[stats::complete.cases(db),] 

for(i in 1:k){
colnames(db)[i]<-colnames(r_t[[i]])
}

TT<-nrow(db)

db<-zoo::coredata(db)

##### matrices and vectors

S<-stats::cov(db)
H_t<-array(S,dim=c(k,k,TT))

##### cycle

for(tt in (V+1):TT){
X<-db[(tt-1):(tt-V),]
arr<-array(apply(X, 1, function (x)  tcrossprod(x)), dim=c(k,k,V))
H_t[,,tt]<-apply( arr , 1:2 , mean )

}

results<-list(
"H_t"=H_t
)

return(results)

}

#' sBEKK log-likelihood 
#'
#' Obtains the log-likelihood of the scalar BEKK model.
#' @param param Vector of starting values. 
#' @param ret Txk matrix of daily returns. At the moment, k can be at most 5
#' @return The resulting vector is the log-likelihood value for each \eqn{t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @keywords internal
#' @export

sBEKK_loglik<-function(param,ret){

##### index

k<-ncol(ret)		# number of assets
TT<-nrow(ret)		# number of daily observations

if(k==2){

c_11<-param[1]
c_21<-param[2]
c_22<-param[3]
a<-param[4]
b<-param[5]

C_mat<-matrix(0,ncol=k,nrow=k)

C_mat[1,1]<-c_11
C_mat[2,1]<-c_21
C_mat[2,2]<-c_22

} else if (k==3){

C_mat<-matrix(0,ncol=k,nrow=k)

c_11<-param[1]
c_21<-param[2]
c_22<-param[3]
c_31<-param[4]
c_32<-param[5]
c_33<-param[6]
a<-param[7]
b<-param[8]

C_mat[1,1]<-c_11
C_mat[2,1]<-c_21
C_mat[2,2]<-c_22
C_mat[3,1]<-c_31
C_mat[3,2]<-c_32
C_mat[3,3]<-c_33

} else if (k==4){

C_mat<-matrix(0,ncol=k,nrow=k)

c_11<-param[1]
c_21<-param[2]
c_22<-param[3]
c_31<-param[4]
c_32<-param[5]
c_33<-param[6]
c_41<-param[7]
c_42<-param[8]
c_43<-param[9]
c_44<-param[10]
a<-param[11]
b<-param[12]

C_mat[1,1]<-c_11
C_mat[2,1]<-c_21
C_mat[2,2]<-c_22
C_mat[3,1]<-c_31
C_mat[3,2]<-c_32
C_mat[3,3]<-c_33
C_mat[4,1]<-c_41
C_mat[4,2]<-c_42
C_mat[4,3]<-c_43
C_mat[4,4]<-c_44

} else if (k==5){

C_mat<-matrix(0,ncol=k,nrow=k)

c_11<-param[1]
c_21<-param[2]
c_22<-param[3]
c_31<-param[4]
c_32<-param[5]
c_33<-param[6]
c_41<-param[7]
c_42<-param[8]
c_43<-param[9]
c_44<-param[10]
c_51<-param[11]
c_52<-param[12]
c_53<-param[13]
c_54<-param[14]
c_55<-param[15]
a<-param[16]
b<-param[17]

C_mat[1,1]<-c_11
C_mat[2,1]<-c_21
C_mat[2,2]<-c_22
C_mat[3,1]<-c_31
C_mat[3,2]<-c_32
C_mat[3,3]<-c_33
C_mat[4,1]<-c_41
C_mat[4,2]<-c_42
C_mat[4,3]<-c_43
C_mat[4,4]<-c_44
C_mat[5,1]<-c_51
C_mat[5,2]<-c_52
C_mat[5,3]<-c_53
C_mat[5,4]<-c_54
C_mat[5,5]<-c_55

}

W<-C_mat%*%t(C_mat)

##### matrices and vectors

H_0<-stats::cov(ret)

H_t<-array(NA,dim=c(k,k,TT))
H_t[,,1]<-H_0
log_det_H_t<-rep(0,TT)

ret_H_t_solved_ret<-rep(0,TT)


################# likelihood

ll<-rep(0,TT)

for(tt in 2:TT){
Cross_prod<-tcrossprod(ret[tt-1,],ret[tt-1,])
H_t[,,tt]<-W+a^2*Cross_prod+b^2*H_t[,,tt-1]
log_det_H_t[tt]<-log(Det(H_t[,,tt]))
H_t_solved<-Inv(H_t[,,tt])
ret_H_t_solved_ret[tt]<-(ret[tt,])%*%H_t_solved%*%cbind(ret[tt,])
}

ll<- - (log_det_H_t+ret_H_t_solved_ret)

#sum(ll)

return(ll)

}

#' sBEKK covariance matrix
#'
#' Obtains the conditional covariance matrix from the scalar BEKK model.
#' @param param Vector of starting values. 
#' @param ret Txk matrix of daily returns. At the moment, k can be at most 5
#' @return The resulting vector is the log-likelihood value for each \eqn{t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @keywords internal
#' @export

sBEKK_mat_est<-function(param,ret){


##### index

k<-ncol(ret)		# number of assets
TT<-nrow(ret)		# number of daily observations

if(k==2){

c_11<-param[1]
c_21<-param[2]
c_22<-param[3]
a<-param[4]
b<-param[5]

C_mat<-matrix(0,ncol=k,nrow=k)

C_mat[1,1]<-c_11
C_mat[2,1]<-c_21
C_mat[2,2]<-c_22

} else if (k==3){

C_mat<-matrix(0,ncol=k,nrow=k)

c_11<-param[1]
c_21<-param[2]
c_22<-param[3]
c_31<-param[4]
c_32<-param[5]
c_33<-param[6]
a<-param[7]
b<-param[8]

C_mat[1,1]<-c_11
C_mat[2,1]<-c_21
C_mat[2,2]<-c_22
C_mat[3,1]<-c_31
C_mat[3,2]<-c_32
C_mat[3,3]<-c_33

} else if (k==4){

C_mat<-matrix(0,ncol=k,nrow=k)

c_11<-param[1]
c_21<-param[2]
c_22<-param[3]
c_31<-param[4]
c_32<-param[5]
c_33<-param[6]
c_41<-param[7]
c_42<-param[8]
c_43<-param[9]
c_44<-param[10]
a<-param[11]
b<-param[12]

C_mat[1,1]<-c_11
C_mat[2,1]<-c_21
C_mat[2,2]<-c_22
C_mat[3,1]<-c_31
C_mat[3,2]<-c_32
C_mat[3,3]<-c_33
C_mat[4,1]<-c_41
C_mat[4,2]<-c_42
C_mat[4,3]<-c_43
C_mat[4,4]<-c_44


} else if (k==5){

C_mat<-matrix(0,ncol=k,nrow=k)

c_11<-param[1]
c_21<-param[2]
c_22<-param[3]
c_31<-param[4]
c_32<-param[5]
c_33<-param[6]
c_41<-param[7]
c_42<-param[8]
c_43<-param[9]
c_44<-param[10]
c_51<-param[11]
c_52<-param[12]
c_53<-param[13]
c_54<-param[14]
c_55<-param[15]
a<-param[16]
b<-param[17]

C_mat[1,1]<-c_11
C_mat[2,1]<-c_21
C_mat[2,2]<-c_22
C_mat[3,1]<-c_31
C_mat[3,2]<-c_32
C_mat[3,3]<-c_33
C_mat[4,1]<-c_41
C_mat[4,2]<-c_42
C_mat[4,3]<-c_43
C_mat[4,4]<-c_44
C_mat[5,1]<-c_51
C_mat[5,2]<-c_52
C_mat[5,3]<-c_53
C_mat[5,4]<-c_54
C_mat[5,5]<-c_55

}


W<-C_mat%*%t(C_mat)

##### matrices and vectors

H_0<-stats::cov(ret)

H_t<-array(NA,dim=c(k,k,TT))
H_t[,,1]<-H_0

################# cov-est

for(tt in 2:TT){
Cross_prod<-tcrossprod(ret[tt-1,],ret[tt-1,])
H_t[,,tt]<-W+a^2*Cross_prod+b^2*H_t[,,tt-1]
}

return(H_t)

}


#' dBEKK log-likelihood 
#'
#' Obtains the log-likelihood of the diagonal BEKK model.
#' @param param Vector of starting values. 
#' @param ret Txk matrix of daily returns. At the moment, k can be at most 5
#' @return The resulting vector is the log-likelihood value for each \eqn{t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @keywords internal
#' @export

dBEKK_loglik<-function(param,ret){

##### index

k<-ncol(ret)		# number of assets
TT<-nrow(ret)		# number of daily observations

if(k==2){

c_11<-param[1]
c_21<-param[2]
c_22<-param[3]
a_11<-param[4]
a_22<-param[5]
b_11<-param[6]
b_22<-param[7]

A_mat<-B_mat<-C_mat<-matrix(0,ncol=k,nrow=k)

C_mat[1,1]<-c_11
C_mat[2,1]<-c_21
C_mat[2,2]<-c_22

A_mat[1,1]<-a_11
A_mat[2,2]<-a_22

B_mat[1,1]<-b_11
B_mat[2,2]<-b_22


} else if (k==3){

A_mat<-B_mat<-C_mat<-matrix(0,ncol=k,nrow=k)

c_11<-param[1]
c_21<-param[2]
c_22<-param[3]
c_31<-param[4]
c_32<-param[5]
c_33<-param[6]
a_11<-param[7]
a_22<-param[8]
a_33<-param[9]
b_11<-param[10]
b_22<-param[11]
b_33<-param[12]

C_mat[1,1]<-c_11
C_mat[2,1]<-c_21
C_mat[2,2]<-c_22
C_mat[3,1]<-c_31
C_mat[3,2]<-c_32
C_mat[3,3]<-c_33


A_mat[1,1]<-a_11
A_mat[2,2]<-a_22
A_mat[3,3]<-a_33

B_mat[1,1]<-b_11
B_mat[2,2]<-b_22
B_mat[3,3]<-b_33

} else if (k==4){

A_mat<-B_mat<-C_mat<-matrix(0,ncol=k,nrow=k)

c_11<-param[1]
c_21<-param[2]
c_22<-param[3]
c_31<-param[4]
c_32<-param[5]
c_33<-param[6]
c_41<-param[7]
c_42<-param[8]
c_43<-param[9]
c_44<-param[10]

a_11<-param[11]
a_22<-param[12]
a_33<-param[13]
a_44<-param[14]

b_11<-param[15]
b_22<-param[16]
b_33<-param[17]
b_44<-param[18]

C_mat[1,1]<-c_11
C_mat[2,1]<-c_21
C_mat[2,2]<-c_22
C_mat[3,1]<-c_31
C_mat[3,2]<-c_32
C_mat[3,3]<-c_33
C_mat[4,1]<-c_41
C_mat[4,2]<-c_42
C_mat[4,3]<-c_43
C_mat[4,4]<-c_44

A_mat[1,1]<-a_11
A_mat[2,2]<-a_22
A_mat[3,3]<-a_33
A_mat[4,4]<-a_44

B_mat[1,1]<-b_11
B_mat[2,2]<-b_22
B_mat[3,3]<-b_33
B_mat[4,4]<-b_44

} else if (k==5){

A_mat<-B_mat<-C_mat<-matrix(0,ncol=k,nrow=k)

c_11<-param[1]
c_21<-param[2]
c_22<-param[3]
c_31<-param[4]
c_32<-param[5]
c_33<-param[6]
c_41<-param[7]
c_42<-param[8]
c_43<-param[9]
c_44<-param[10]
c_51<-param[11]
c_52<-param[12]
c_53<-param[13]
c_54<-param[14]
c_55<-param[15]

a_11<-param[16]
a_22<-param[17]
a_33<-param[18]
a_44<-param[19]
a_55<-param[20]

b_11<-param[21]
b_22<-param[22]
b_33<-param[23]
b_44<-param[24]
b_55<-param[25]

C_mat[1,1]<-c_11
C_mat[2,1]<-c_21
C_mat[2,2]<-c_22
C_mat[3,1]<-c_31
C_mat[3,2]<-c_32
C_mat[3,3]<-c_33
C_mat[4,1]<-c_41
C_mat[4,2]<-c_42
C_mat[4,3]<-c_43
C_mat[4,4]<-c_44
C_mat[5,1]<-c_51
C_mat[5,2]<-c_52
C_mat[5,3]<-c_53
C_mat[5,4]<-c_54
C_mat[5,5]<-c_55

A_mat[1,1]<-a_11
A_mat[2,2]<-a_22
A_mat[3,3]<-a_33
A_mat[4,4]<-a_44
A_mat[5,5]<-a_55

B_mat[1,1]<-b_11
B_mat[2,2]<-b_22
B_mat[3,3]<-b_33
B_mat[4,4]<-b_44
B_mat[5,5]<-b_55

}


W<-C_mat%*%t(C_mat)

##### matrices and vectors

H_0<-stats::cov(ret)

H_t<-array(NA,dim=c(k,k,TT))
H_t[,,1]<-H_0
log_det_H_t<-rep(0,TT)

ret_H_t_solved_ret<-rep(0,TT)


################# likelihood

ll<-rep(0,TT)

for(tt in 2:TT){
Cross_prod<-tcrossprod(ret[tt-1,],ret[tt-1,])
H_t[,,tt]<-W+A_mat%*%Cross_prod%*%t(A_mat)+B_mat%*%H_t[,,tt-1]%*%t(B_mat)
log_det_H_t[tt]<-log(Det(H_t[,,tt]))
H_t_solved_1<-suppressWarnings(
tryCatch(Inv(H_t[,,tt]), error = function(e) return(NA))
)
if (all(is.na(H_t_solved_1))){
ret_H_t_solved_ret[tt] <- Inf
} else{
H_t_solved<-H_t_solved_1
ret_H_t_solved_ret[tt]<-(ret[tt,])%*%H_t_solved%*%cbind(ret[tt,])
}

}


ll<- - (log_det_H_t+ret_H_t_solved_ret)

#sum(ll)

return(ll)

}

#' dBEKK covariance matrix
#'
#' Obtains the conditional covariance matrix from the diagonal BEKK model.
#' @param param Vector of starting values. 
#' @param ret Txk matrix of daily returns. At the moment, k can be at most 4
#' @return The resulting vector is the log-likelihood value for each \eqn{t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @keywords internal
#' @export

dBEKK_mat_est<-function(param,ret){


##### index

k<-ncol(ret)		# number of assets
TT<-nrow(ret)		# number of daily observations

if(k==2){

c_11<-param[1]
c_21<-param[2]
c_22<-param[3]
a_11<-param[4]
a_22<-param[5]
b_11<-param[6]
b_22<-param[7]

A_mat<-B_mat<-C_mat<-matrix(0,ncol=k,nrow=k)

C_mat[1,1]<-c_11
C_mat[2,1]<-c_21
C_mat[2,2]<-c_22

A_mat[1,1]<-a_11
A_mat[2,2]<-a_22

B_mat[1,1]<-b_11
B_mat[2,2]<-b_22


} else if (k==3){

A_mat<-B_mat<-C_mat<-matrix(0,ncol=k,nrow=k)

c_11<-param[1]
c_21<-param[2]
c_22<-param[3]
c_31<-param[4]
c_32<-param[5]
c_33<-param[6]
a_11<-param[7]
a_22<-param[8]
a_33<-param[9]
b_11<-param[10]
b_22<-param[11]
b_33<-param[12]

C_mat[1,1]<-c_11
C_mat[2,1]<-c_21
C_mat[2,2]<-c_22
C_mat[3,1]<-c_31
C_mat[3,2]<-c_32
C_mat[3,3]<-c_33


A_mat[1,1]<-a_11
A_mat[2,2]<-a_22
A_mat[3,3]<-a_33

B_mat[1,1]<-b_11
B_mat[2,2]<-b_22
B_mat[3,3]<-b_33

} else if (k==4){

A_mat<-B_mat<-C_mat<-matrix(0,ncol=k,nrow=k)

c_11<-param[1]
c_21<-param[2]
c_22<-param[3]
c_31<-param[4]
c_32<-param[5]
c_33<-param[6]
c_41<-param[7]
c_42<-param[8]
c_43<-param[9]
c_44<-param[10]

a_11<-param[11]
a_22<-param[12]
a_33<-param[13]
a_44<-param[14]

b_11<-param[15]
b_22<-param[16]
b_33<-param[17]
b_44<-param[18]

C_mat[1,1]<-c_11
C_mat[2,1]<-c_21
C_mat[2,2]<-c_22
C_mat[3,1]<-c_31
C_mat[3,2]<-c_32
C_mat[3,3]<-c_33
C_mat[4,1]<-c_41
C_mat[4,2]<-c_42
C_mat[4,3]<-c_43
C_mat[4,4]<-c_44

A_mat[1,1]<-a_11
A_mat[2,2]<-a_22
A_mat[3,3]<-a_33
A_mat[4,4]<-a_44

B_mat[1,1]<-b_11
B_mat[2,2]<-b_22
B_mat[3,3]<-b_33
B_mat[4,4]<-b_44

} else if (k==5){

A_mat<-B_mat<-C_mat<-matrix(0,ncol=k,nrow=k)

c_11<-param[1]
c_21<-param[2]
c_22<-param[3]
c_31<-param[4]
c_32<-param[5]
c_33<-param[6]
c_41<-param[7]
c_42<-param[8]
c_43<-param[9]
c_44<-param[10]
c_51<-param[11]
c_52<-param[12]
c_53<-param[13]
c_54<-param[14]
c_55<-param[15]

a_11<-param[16]
a_22<-param[17]
a_33<-param[18]
a_44<-param[19]
a_55<-param[20]

b_11<-param[21]
b_22<-param[22]
b_33<-param[23]
b_44<-param[24]
b_55<-param[25]

C_mat[1,1]<-c_11
C_mat[2,1]<-c_21
C_mat[2,2]<-c_22
C_mat[3,1]<-c_31
C_mat[3,2]<-c_32
C_mat[3,3]<-c_33
C_mat[4,1]<-c_41
C_mat[4,2]<-c_42
C_mat[4,3]<-c_43
C_mat[4,4]<-c_44
C_mat[5,1]<-c_51
C_mat[5,2]<-c_52
C_mat[5,3]<-c_53
C_mat[5,4]<-c_54
C_mat[5,5]<-c_55

A_mat[1,1]<-a_11
A_mat[2,2]<-a_22
A_mat[3,3]<-a_33
A_mat[4,4]<-a_44
A_mat[5,5]<-a_55

B_mat[1,1]<-b_11
B_mat[2,2]<-b_22
B_mat[3,3]<-b_33
B_mat[4,4]<-b_44
B_mat[5,5]<-b_55

}


W<-C_mat%*%t(C_mat)

##### matrices and vectors

H_0<-stats::cov(ret)

H_t<-array(NA,dim=c(k,k,TT))
H_t[,,1]<-H_0

################# cov-est

for(tt in 2:TT){
Cross_prod<-tcrossprod(ret[tt-1,],ret[tt-1,])
H_t[,,tt]<-W+A_mat%*%Cross_prod%*%t(A_mat)+B_mat%*%H_t[,,tt-1]%*%t(B_mat)
}

return(H_t)

}


#' cDCC log-likelihood (second step)
#'
#' Obtains the log-likelihood of the cDCC model in the second step.
#' For details, see \insertCite{aielli2013dynamic;textual}{dccmidas} and \insertCite{dcc_engle_2002;textual}{dccmidas}.
#' @param param Vector of starting values. 
#' @param res Array of standardized daily returns, coming from the first step estimation.
#' @param K_c **optional** Number of initial observations to exclude from the estimation
#' @return The resulting vector is the log-likelihood value for each \eqn{t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @import rumidas
#' @references
#' \insertAllCited{} 
#' @keywords internal
#' @export

dcc_loglik<-function(param,res,K_c=NULL){

a<-param[1]
b<-param[2]

##### index

Num_col<-dim(res)[1]		# number of assets
TT<-dim(res)[3]			# number of daily observations

##### K_c

if(is.null(K_c)){
K_c<-2
} else {
K_c<-K_c+1
}

##### matrices and vectors

Prod_eps_t<-array(0,dim=c(Num_col,Num_col,TT))
Q_t<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
Q_t_star<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
R_t<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
log_det_R_t<-rep(0,TT)
R_t_solved<-array(1,dim=c(Num_col,Num_col,TT))
Eps_t_cross_prod<-rep(0,TT)
Eps_t_R_t_Eps_t<-rep(0,TT)


################# likelihood

ll<-rep(0,TT)

S<-stats::cov(t(apply(res, 3L, c)))

for(tt in K_c:TT){
Q_t_star<-sqrt(diag(diag(Q_t[,,tt-1])))
Q_t[,,tt]<-(1-a-b)*S + a*(Q_t_star%*%res[,,tt-1]%*%t(res[,,tt-1])%*%Q_t_star) + b*Q_t[,,tt-1]
Q_t_star_inv<-Inv(sqrt(diag(diag(Q_t[,,tt]))))
R_t[,,tt]<-Q_t_star_inv%*%Q_t[,,tt]%*%Q_t_star_inv
log_det_R_t[tt]<-log(Det(R_t[,,tt]))
R_t_solved[,,tt]<-Inv(R_t[,,tt])
#Eps_t_cross_prod[tt]<-rbind(res[,,tt])%*%cbind(res[,,tt])
Eps_t_R_t_Eps_t[tt]<-rbind(res[,,tt])%*%R_t_solved[,,tt]%*%cbind(res[,,tt])
}

ll<- - (log_det_R_t+Eps_t_R_t_Eps_t)

#sum(ll)

return(ll)

}

#' Obtains the matrix H_t and R_t, under the cDCC model
#'
#' Obtains the matrix H_t and R_t, under the cDCC model
#' For details, see \insertCite{aielli2013dynamic;textual}{dccmidas} and \insertCite{dcc_engle_2002;textual}{dccmidas}.
#' @param est_param Vector of estimated values 
#' @param res Array of standardized daily returns, coming from the first step estimation
#' @param Dt Diagonal matrix of standard deviations
#' @param K_c **optional** Number of initial observations to exclude from the H_t and R_t calculation
#' @return A list with the \eqn{H_t} and \eqn{R_t} matrices, for each \eqn{t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @keywords internal
#' @export

dcc_mat_est<-function(est_param,res,Dt,K_c){

a<-est_param[1]
b<-est_param[2]

##### index

Num_col<-dim(res)[1]		# number of assets
TT<-dim(res)[3]			# number of daily observations

##### K_c

if(is.null(K_c)){
K_c<-2
} else {
K_c<-K_c+1
}


##### matrices and vectors

Prod_eps_t<-array(0,dim=c(Num_col,Num_col,TT))
Q_t<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
Q_t_star<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
Q_t_star_inv<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))

R_t<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))

S<-stats::cov(t(apply(res, 3L, c)))
H_t<-array(S,dim=c(Num_col,Num_col,TT))

################# H_t and R_t


for(tt in K_c:TT){
Q_t_star<-sqrt(diag(diag(Q_t[,,tt-1])))
Q_t[,,tt]<-(1-a-b)*S + a*(Q_t_star%*%res[,,tt-1]%*%t(res[,,tt-1])%*%Q_t_star) + b*Q_t[,,tt-1]
Q_t_star_inv<-Inv(sqrt(diag(diag(Q_t[,,tt]))))
R_t[,,tt]<-Q_t_star_inv%*%Q_t[,,tt]%*%Q_t_star_inv
H_t[,,tt]<-Dt[,,tt]%*%R_t[,,tt]%*%Dt[,,tt]

}

#H_t_m<-apply(H_t[,,2:TT],c(1,2),mean)
#H_t[,,1]<-H_t_m

results<-list(
"H_t"=H_t,
"R_t"=R_t
)

}

#' A-DCC log-likelihood (second step)
#'
#' Obtains the log-likelihood of the A-DCC model in the second step.
#' For details, see \insertCite{cappiello2006asymmetric;textual}{dccmidas} and \insertCite{dcc_engle_2002;textual}{dccmidas}.
#' @param param Vector of starting values. 
#' @param res Array of standardized daily returns, coming from the first step estimation.
#' @param K_c **optional** Number of initial observations to exclude from the estimation
#' @return The resulting vector is the log-likelihood value for each \eqn{t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @keywords internal
#' @export

a_dcc_loglik<-function(param,res,K_c=NULL){

a<-param[1]
b<-param[2]
g<-param[3]

##### index

Num_col<-dim(res)[1]		# number of assets
TT<-dim(res)[3]			# number of daily observations

##### K_c

if(is.null(K_c)){
K_c<-2
} else {
K_c<-K_c+1
}


##### matrices and vectors

Prod_eps_t<-array(0,dim=c(Num_col,Num_col,TT))
Q_t<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
Q_t_star<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
R_t<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
log_det_R_t<-rep(0,TT)
R_t_solved<-array(1,dim=c(Num_col,Num_col,TT))
Eps_t_cross_prod<-rep(0,TT)
Eps_t_R_t_Eps_t<-rep(0,TT)

eta<-array(0,dim=c(Num_col,1,TT))
eta<-ifelse(res<0,1,0)*res

################# likelihood

ll<-rep(0,TT)

S<-stats::cov(t(apply(res, 3L, c)))
N<-stats::cov(t(apply(eta, 3L, c)))

for(tt in (K_c):TT){

Q_t[,,tt]<-(1-a-b)*S - g*N + a*(cbind(res[,,tt-1])%*%rbind(res[,,tt-1])) + 
b*Q_t[,,tt-1] + g*eta[,,tt-1]%*%t(eta[,,tt-1])
Q_t_star_inv<-Inv(sqrt(diag(diag(Q_t[,,tt]))))
R_t[,,tt]<-Q_t_star_inv%*%Q_t[,,tt]%*%Q_t_star_inv
log_det_R_t[tt]<-log(Det(R_t[,,tt]))
R_t_solved[,,tt]<-Inv(R_t[,,tt])
Eps_t_R_t_Eps_t[tt]<-rbind(res[,,tt])%*%R_t_solved[,,tt]%*%cbind(res[,,tt])
}

ll<- - (log_det_R_t+Eps_t_R_t_Eps_t)

return(ll)

}

#' Obtains the matrix H_t and R_t, under the A-DCC model
#'
#' Obtains the matrix H_t and R_t, under the A-DCC model
#' For details, see \insertCite{cappiello2006asymmetric;textual}{dccmidas} and \insertCite{dcc_engle_2002;textual}{dccmidas}.
#' @param est_param Vector of estimated values
#' @param res Array of standardized daily returns, coming from the first step estimation
#' @param Dt Diagonal matrix of standard deviations
#' @param K_c **optional** Number of initial observations to exclude from the H_t and R_t calculation
#' @return A list with the \eqn{H_t} and \eqn{R_t} matrices, for each \eqn{t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @keywords internal
#' @export

a_dcc_mat_est<-function(est_param,res,Dt,K_c=NULL){

a<-est_param[1]
b<-est_param[2]
g<-est_param[3]

##### index

Num_col<-dim(res)[1]		# number of assets
TT<-dim(res)[3]			# number of daily observations

##### K_c

if(is.null(K_c)){
K_c<-2
} else {
K_c<-K_c+1
}

##### matrices and vectors

Prod_eps_t<-array(0,dim=c(Num_col,Num_col,TT))
Q_t<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
Q_t_star<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
Q_t_star_inv<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))

R_t<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))

S<-stats::cov(t(apply(res, 3L, c)))
H_t<-array(S,dim=c(Num_col,Num_col,TT))

eta<-array(0,dim=c(Num_col,1,TT))
eta<-ifelse(res<0,1,0)*res

################# H_t and R_t

S<-stats::cov(t(apply(res, 3L, c)))
N<-stats::cov(t(apply(eta, 3L, c)))

for(tt in (K_c):TT){
Q_t[,,tt]<-(1-a-b)*S - g*N + a*(cbind(res[,,tt-1])%*%rbind(res[,,tt-1])) + 
b*Q_t[,,tt-1] + g*eta[,,tt-1]%*%t(eta[,,tt-1])
Q_t_star_inv<-Inv(sqrt(diag(diag(Q_t[,,tt]))))
R_t[,,tt]<-Q_t_star_inv%*%Q_t[,,tt]%*%Q_t_star_inv
H_t[,,tt]<-Dt[,,tt]%*%R_t[,,tt]%*%Dt[,,tt]
}

#H_t_m<-apply(H_t[,,(2):TT],c(1,2),mean)
#H_t[,,1]<-H_t_m

results<-list(
"H_t"=H_t,
"R_t"=R_t
)

}

#' A-DCC-MIDAS log-likelihood (second step)
#'
#' Obtains the log-likelihood of the A-DCC-MIDAS model in the second step.
#' For details, see \insertCite{cappiello2006asymmetric;textual}{dccmidas} and \insertCite{dcc_engle_2002;textual}{dccmidas}.
#' @param param Vector of starting values. 
#' @param res Array of standardized daily returns, coming from the first step estimation.
#' @return The resulting vector is the log-likelihood value for each \eqn{t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @import rumidas
#' @references
#' \insertAllCited{} 
#' @keywords internal
#' @export

a_dccmidas_loglik<-function(param,res,lag_fun="Beta",N_c,K_c){

a<-param[1]
b<-param[2]
g<-param[3]
w1<- ifelse(lag_fun=="Beta",1,0)
w2<-param[4]

##### index

Num_col<-dim(res)[1]		# number of assets
TT<-dim(res)[3]			# number of daily observations

##### matrices and vectors

C_t<-array(0,dim=c(Num_col,Num_col,TT))
V_t<-array(0,dim=c(Num_col,Num_col,TT))
Prod_eps_t<-array(0,dim=c(Num_col,Num_col,TT))
Q_t<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
Q_t_star<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
R_t<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
log_det_R_t<-rep(0,TT)
R_t_solved<-array(1,dim=c(Num_col,Num_col,TT))
Eps_t_cross_prod<-rep(0,TT)
Eps_t_R_t_Eps_t<-rep(0,TT)

eta<-array(0,dim=c(Num_col,1,TT))
eta<-ifelse(res<0,1,0)*res

##### first cycle

for(tt in (N_c+1):TT){
V_t[,,tt]<-rowSums(res[,1,tt:(tt-N_c)]*(res[,1,tt:(tt-N_c)]))*diag(Num_col)
Prod_eps_t[,,tt]<-res[,,tt:(tt-N_c)]%*%t(res[,,tt:(tt-N_c)])
V_t_0.5<-Inv(sqrt(V_t[,,tt]))
C_t[,,tt]<-V_t_0.5%*%Prod_eps_t[,,tt]%*%V_t_0.5
}

#### R_t_bar

weight_fun<-ifelse(lag_fun=="Beta",rumidas::beta_function,rumidas::exp_almon)

betas<-c(rev(weight_fun(1:(K_c+1),(K_c+1),w1,w2))[2:(K_c+1)],0)
R_t_bar<-array(1,dim=c(Num_col,Num_col,TT))

matrix_id<-matrix(1:Num_col^2,ncol=Num_col)
matrix_id_2<-which(matrix_id==1:Num_col^2, arr.ind=TRUE)

for(i in 1:nrow(matrix_id_2)){
R_t_bar[matrix_id_2[i,1],matrix_id_2[i,2],]                      <- suppressWarnings(
roll::roll_sum(C_t[matrix_id_2[i,1],matrix_id_2[i,2],], c(K_c+1),weights = betas)
) 
}

R_t_bar[,,1:K_c]<-diag(Num_col)


################# likelihood

ll<-rep(0,TT)

N<-stats::cov(t(apply(eta, 3L, c)))

for(tt in (2):TT){
Q_t[,,tt]<-(1-a-b)*R_t_bar[,,tt] - g*N + a*(cbind(res[,,tt-1])%*%rbind(res[,,tt-1])) + 
b*Q_t[,,tt-1] + g*eta[,,tt-1]%*%t(eta[,,tt-1])
Q_t_star_inv<-Inv(sqrt(diag(diag(Q_t[,,tt]))))
R_t[,,tt]<-Q_t_star_inv%*%Q_t[,,tt]%*%Q_t_star_inv
log_det_R_t[tt]<-log(Det(R_t[,,tt]))
R_t_solved[,,tt]<-Inv(R_t[,,tt])
Eps_t_R_t_Eps_t[tt]<-rbind(res[,,tt])%*%R_t_solved[,,tt]%*%cbind(res[,,tt])
}

ll<- - (log_det_R_t+Eps_t_R_t_Eps_t)

return(ll)

}

#' Obtains the matrix H_t, R_t and long-run correlations, under the A-DCC-MIDAS model
#'
#' Obtains the matrix H_t, R_t and long-run correlations, under the A-DCC-MIDAS model
#' For details, see \insertCite{colacito2011component;textual}{dccmidas} and \insertCite{dcc_engle_2002;textual}{dccmidas}.
#' @param est_param Vector of estimated values
#' @param res Array of standardized daily returns, coming from the first step estimation
#' @param Dt Matrix of conditional standard deviations (coming from the first step)
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively
#' @param N_c Number of (lagged) realizations to use for the standarized residuals forming the long-run correlation
#' @param K_c Number of (lagged) realizations to use for the long-run correlation
#' @return A list with the \eqn{H_t}, \eqn{R_t} and long-run correlaton matrices, for each \eqn{t}.
#' @importFrom Rdpack reprompt
#' @import roll
#' @import rumidas
#' @references
#' \insertAllCited{} 
#' @keywords internal
#' @export


a_dccmidas_mat_est<-function(est_param,res,Dt,lag_fun="Beta",N_c,K_c){

a<-est_param[1]
b<-est_param[2]
g<-est_param[3]
w1<- ifelse(lag_fun=="Beta",1,0)
w2<-est_param[4]

##### index

Num_col<-dim(res)[1]		# number of assets
TT<-dim(res)[3]			# number of daily observations

##### matrices and vectors

C_t<-array(0,dim=c(Num_col,Num_col,TT))
V_t<-array(0,dim=c(Num_col,Num_col,TT))
Prod_eps_t<-array(0,dim=c(Num_col,Num_col,TT))
Q_t<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
R_t<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
log_det_R_t<-rep(0,TT)
R_t_solved<-array(1,dim=c(Num_col,Num_col,TT))
Eps_t_cross_prod<-rep(0,TT)
Eps_t_R_t_Eps_t<-rep(0,TT)

S<-stats::cov(t(apply(res, 3L, c)))
H_t<-array(S,dim=c(Num_col,Num_col,TT))

##### first cycle

for(tt in (N_c+1):TT){
V_t[,,tt]<-rowSums(res[,1,tt:(tt-N_c)]*(res[,1,tt:(tt-N_c)]))*diag(Num_col)
Prod_eps_t[,,tt]<-res[,,tt:(tt-N_c)]%*%t(res[,,tt:(tt-N_c)])
V_t_0.5<-Inv(sqrt(V_t[,,tt]))
C_t[,,tt]<-V_t_0.5%*%Prod_eps_t[,,tt]%*%V_t_0.5
}

#### R_t_bar

weight_fun<-ifelse(lag_fun=="Beta",beta_function,exp_almon)

betas<-c(rev(weight_fun(1:(K_c+1),(K_c+1),w1,w2))[2:(K_c+1)],0)
R_t_bar<-array(1,dim=c(Num_col,Num_col,TT))

matrix_id<-matrix(1:Num_col^2,ncol=Num_col)
matrix_id_2<-which(matrix_id==1:Num_col^2, arr.ind=TRUE)

for(i in 1:nrow(matrix_id_2)){
R_t_bar[matrix_id_2[i,1],matrix_id_2[i,2],]                      <- suppressWarnings(
roll_sum(C_t[matrix_id_2[i,1],matrix_id_2[i,2],], c(K_c+1),weights = betas)
) 
}

R_t_bar[,,1:K_c]<-diag(rep(1,Num_col))

################# H_t, R_t and R_t_bar

eta<-array(0,dim=c(Num_col,1,TT))
eta<-ifelse(res<0,1,0)*res


N<-stats::cov(t(apply(eta, 3L, c)))

for(tt in (K_c+1):TT){
Q_t[,,tt]<-(1-a-b)*R_t_bar[,,tt] - g*N + a*(cbind(res[,,tt-1])%*%rbind(res[,,tt-1])) + 
b*Q_t[,,tt-1] + g*eta[,,tt-1]%*%t(eta[,,tt-1])
Q_t_star<-Inv(sqrt(diag(diag(Q_t[,,tt]))))
R_t[,,tt]<-Q_t_star%*%Q_t[,,tt]%*%Q_t_star
H_t[,,tt]<-Dt[,,tt]%*%R_t[,,tt]%*%Dt[,,tt]
}

H_t_m<-apply(H_t[,,(K_c+1):TT],c(1,2),mean)
H_t[,,1:K_c]<-H_t_m

results<-list(
"H_t"=H_t,
"R_t"=R_t,
"R_t_bar"=R_t_bar
)

return(results)

}

#' DECO log-likelihood (second step)
#'
#' Obtains the log-likelihood of the DECO models in the second step.
#' For details, see \insertCite{engle2012dynamic;textual}{dccmidas}.
#' @param param Vector of starting values. 
#' @param res Array of standardized daily returns, coming from the first step estimation.
#' @param K_c **optional** Number of initial observations to exclude from the estimation
#' @return The resulting vector is the log-likelihood value for each \eqn{t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @keywords internal
#' @export

deco_loglik<-function(param,res,K_c=NULL){

a<-param[1]
b<-param[2]

##### index

Num_col<-dim(res)[1]		# number of assets
TT<-dim(res)[3]			# number of daily observations

##### K_c

if(is.null(K_c)){
K_c<-2
} else {
K_c<-K_c+1
}

##### matrices and vectors

Prod_eps_t<-array(0,dim=c(Num_col,Num_col,TT))
Q_t<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
Q_t_star<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
I_t<-diag(rep(1,Num_col))
J_t<-matrix(rep(1,Num_col^2),ncol=Num_col)
R_t<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
R_t_deco<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
rho_t_deco<-rep(0,TT)
log_det_R_t<-rep(0,TT)
R_t_solved<-array(1,dim=c(Num_col,Num_col,TT))
Eps_t_cross_prod<-rep(0,TT)
Eps_t_R_t_Eps_t<-rep(0,TT)

first_elem<-(Num_col*(Num_col-1))^(-1)
iota<-cbind(rep(1,Num_col))

################# likelihood

ll<-rep(0,TT)

S<-stats::cov(t(apply(res, 3L, c)))

for(tt in K_c:TT){
Q_t_star<-sqrt(diag(diag(Q_t[,,tt-1])))
Q_t[,,tt]<-(1-a-b)*S + a*(Q_t_star%*%res[,,tt-1]%*%t(res[,,tt-1])%*%Q_t_star) + b*Q_t[,,tt-1]
Q_t_star_inv<-Inv(sqrt(diag(diag(Q_t[,,tt]))))
R_t[,,tt]<-Q_t_star_inv%*%Q_t[,,tt]%*%Q_t_star_inv

rho_t_deco<-as.numeric(first_elem*(t(iota)%*%R_t[,,tt]%*%iota - Num_col))

R_t_deco[,,tt]<-(1-rho_t_deco)*I_t+rho_t_deco*J_t

log_det_R_t[tt]<-log(Det(R_t_deco[,,tt]))
R_t_solved[,,tt]<-Inv(R_t_deco[,,tt])

#Eps_t_cross_prod[tt]<-rbind(res[,,tt])%*%cbind(res[,,tt])

Eps_t_R_t_Eps_t[tt]<-rbind(res[,,tt])%*%R_t_solved[,,tt]%*%cbind(res[,,tt])

}

ll<- - (log_det_R_t+Eps_t_R_t_Eps_t)

#sum(ll)

return(ll)

}

#' Obtains the matrix H_t and R_t, under the DECO model
#'
#' Obtains the matrix H_t and R_t, under the DECO model
#' For details, see \insertCite{engle2012dynamic;textual}{dccmidas}.
#' @param est_param Vector of estimated values
#' @param res Array of standardized daily returns, coming from the first step estimation
#' @param Dt Diagonal matrix of standard deviations
#' @param K_c **optional** Number of initial observations to exclude from the H_t and R_t calculation
#' @return A list with the \eqn{H_t} and \eqn{R_t} matrices, for each \eqn{t}.
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @keywords internal
#' @export

deco_mat_est<-function(est_param,res,Dt,K_c=NULL){

a<-est_param[1]
b<-est_param[2]

##### index

Num_col<-dim(res)[1]		# number of assets
TT<-dim(res)[3]			# number of daily observations

##### K_c

if(is.null(K_c)){
K_c<-2
} else {
K_c<-K_c+1
}

##### matrices and vectors

Prod_eps_t<-array(0,dim=c(Num_col,Num_col,TT))
Q_t<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
Q_t_star<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
I_t<-diag(rep(1,Num_col))
J_t<-matrix(rep(1,Num_col^2),ncol=Num_col)
R_t<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
R_t_deco<-array(diag(rep(1,Num_col)),dim=c(Num_col,Num_col,TT))
rho_t_deco<-rep(0,TT)
log_det_R_t<-rep(0,TT)
R_t_solved<-array(1,dim=c(Num_col,Num_col,TT))
Eps_t_cross_prod<-rep(0,TT)
Eps_t_R_t_Eps_t<-rep(0,TT)

first_elem<-(Num_col*(Num_col-1))^(-1)
iota<-cbind(rep(1,Num_col))

################# H_t and R_t

S<-stats::cov(t(apply(res, 3L, c)))
H_t<-array(S,dim=c(Num_col,Num_col,TT))

for(tt in (K_c):TT){
Q_t_star<-sqrt(diag(diag(Q_t[,,tt-1])))
Q_t[,,tt]<-(1-a-b)*S + a*(Q_t_star%*%res[,,tt-1]%*%t(res[,,tt-1])%*%Q_t_star) + b*Q_t[,,tt-1]
Q_t_star_inv<-Inv(sqrt(diag(diag(Q_t[,,tt]))))
R_t[,,tt]<-Q_t_star_inv%*%Q_t[,,tt]%*%Q_t_star_inv

rho_t_deco<-as.numeric(first_elem*(t(iota)%*%R_t[,,tt]%*%iota - Num_col))

R_t_deco[,,tt]<-(1-rho_t_deco)*I_t+rho_t_deco*J_t

H_t[,,tt]<-Dt[,,tt]%*%R_t_deco[,,tt]%*%Dt[,,tt]

}

H_t_m<-apply(H_t[,,(2):TT],c(1,2),mean)
H_t[,,1]<-H_t_m


results<-list(
"H_t"=H_t,
"R_t"=R_t_deco
)

}

#' Standard errors for the Quasi Maximum Likelihood estimator
#'
#' Obtains the standard errors for the Quasi Maximum Likelihood (QML) estimator.
#' @param est It is the output of the maximum likelihood estimation process.
#' @return The resulting vector represents the QML standard errors.
#' @importFrom Rdpack reprompt
#' @import maxLik
#' @keywords internal

QMLE_sd<-function(est){

############### QMLE standard errors (-H)^-1 OPG (-H)^-1

H_hat_inv<- solve(-maxLik::hessian(est))
OPG<-t(est$gradientObs)%*%est$gradientObs

return(diag(H_hat_inv%*%OPG%*%H_hat_inv)^0.5)
}

#' Print method for 'dccmidas' class
#'
#' @param x An object of class 'dccmidas'.
#' @param ... Further arguments passed to or from other methods.
#' @keywords internal
#' @return No return value, called for side effects
#' @export print.dccmidas 
#' @export

print.dccmidas <- function(x, ...) {

options(scipen = 999)

asset_names<-x$assets
model<-x$model
univ_mat<-x$est_univ_model
mat_coef<-x$corr_coef_mat
mult_model<-x$mult_model
est_time<-x$est_time

est_time_print<-paste(
"Estimation time: ",
round(as.numeric(est_time),3),
attributes(est_time)$units)



############################ changing names

if(model=="GM_noskew"){
model="GARCH-MIDAS" 
} else if (model=="GM_skew"){
model="GJR-GARCH-MIDAS" 
} else if (model=="DAGM_skew"){
model="DAGM" 
}  else if (model=="DAGM_noskew"){
model="Asymmetric-GARCH-MIDAS" 
} 

if(mult_model=="DCCMIDAS"){
mult_model<-"DCC-MIDAS"
}

###########################


# univariate matrix

univ_mat_coef<-lapply(univ_mat, "[", 1:3)

names(univ_mat_coef) <- asset_names

# correlation coefficients
coef_char<-as.character(round(mat_coef[,1],4))

row_names<-gsub("\\s", " ", format(rownames(mat_coef), width=9))
coef_val<-gsub("\\s", " ", format(coef_char,width=9))

cat(
cat("\n"),
cat("Assets:",paste(asset_names[1:(length(asset_names)-1)],",",sep=""),
asset_names[length(asset_names)],"\n"),
cat("\n"),
cat(paste("Univariate model:",model),"\n"),
cat("\n"),
cat(paste("Correlation model:",mult_model),"\n"),
cat("\n"),
cat(paste("Coefficients of the correlation model: \n",sep="\n")),
cat(row_names, sep=" ", "\n"),
cat(coef_val, sep=" ", "\n"),
cat("\n"),
cat(est_time_print,"\n"),
cat("\n"))
}

#' Summary method for 'dccmidas' class
#'
#' @param object An object of class 'dccmidas', that is the result of a call to \code{\link{dcc_fit}}.
#' @param ... Additional arguments affecting the summary produced.
#' @importFrom utils capture.output
#' @keywords internal
#' @return Returns the printed tables of the univariate and multivariate steps as well as 
#' some additional information about the number of observations, sample period, and information criteria
#' @export summary.dccmidas 
#' @export

summary.dccmidas <- function(object, ...) {

model<-object$model

####################################### changing names

if(model=="GM_noskew"){
model<-"GARCH-MIDAS" 
} else if (model=="GM_skew"){
model<-"GJR-GARCH-MIDAS" 
} else if (model=="DAGM_skew"){
model<-"DAGM" 
}  else if (model=="DAGM_noskew"){
model<-"Asymmetric-GARCH-MIDAS" 
} 

mult_model<-object$mult_model

if(mult_model=="DCCMIDAS"){
mult_model<-"DCC-MIDAS"
} else if (mult_model=="ADCCMIDAS"){
mult_model<-"A-DCC-MIDAS"
}
#############################################

mat_coef<-object$corr_coef_mat
Obs<-object$obs
Period<-object$period

Period<-paste(substr(Period[1], 1, 10),"/",
substr(Period[2], 1, 10),sep="")

################################ multivariate matrix

p_value<-mat_coef[,4]

sig<-ifelse(p_value<=0.01,"***",ifelse(p_value>0.01&p_value<=0.05,"**",
ifelse(p_value>0.05&p_value<=0.1,"*"," ")))

mat_coef<-round(mat_coef,4)

mat_coef<-cbind(mat_coef,Sig.=sig)

################################ univariate matrix

assets<-object$assets

est_univ<-object$est_univ_model

mat_f<-list()
for(i in 1:length(assets)){
mat<-est_univ[[i]]
p_value<-mat[,4]
sig<-ifelse(p_value<=0.01,"***",ifelse(p_value>0.01&p_value<=0.05,"**",
ifelse(p_value>0.05&p_value<=0.1,"*"," ")))
mat<-round(mat,4)
mat<-data.frame(mat,Sig.=sig)
colnames(mat)[1:4]<-colnames(est_univ[[i]])
mat_f[[i]]<-mat
}
names(mat_f) <- assets

################################# information criteria

llk<-round(object$llk,3)

AIC<-2*nrow(mat_coef) - 2*llk
BIC<-nrow(mat_coef)*log(Obs) - 2*llk

AIC<-round(AIC,3)
BIC<-round(BIC,3)

################################# final print

cat(
cat("\n"),
cat(paste("Univariate model:",model),"\n"),
cat("\n"),
cat("Est. coefficients of the univariate models:\n"),
cat("\n"),
cat(utils::capture.output(mat_f),  sep = '\n'),
cat("---------------------------------------------------","\n"),
cat(paste("Correlation model:",mult_model),"\n"),
cat("\n"),
cat("Est. coefficients of the correlation model:\n"),
cat("\n"),
cat(utils::capture.output(mat_coef),  sep = '\n'),
cat("--- \n"),
cat("Signif. codes: 0.01 '***', 0.05 '**', 0.1 '*' \n"),
cat("\n"),
cat("Obs.:", paste(Obs, ".",sep=""), "Sample Period:", Period, "\n"),
cat("\n"),
cat("LogLik:", paste(llk, ".",sep=""), 
paste("AIC: ", AIC,sep=""), 
paste("BIC: ", BIC,sep=""),
"\n"),
cat("\n"))

}

#' DCC fit (first and second steps)
#'
#' Obtains the estimation of a variety of DCC models, using as univariate models both GARCH and GARCH-MIDAS specifications.
#' @param r_t List of daily returns on which estimate a DCC model. Each daily return must be an 'xts' object.
#' Note that the sample period of the returns should be the same. Otherwise, a merge is performed 
#' @param univ_model Specification of the univariate model. Valid choices are: some of the specifications used in the \code{rugarch} (\link[rugarch]{ugarchspec})
#' and \code{rumidas} (\link[rumidas]{ugmfit}) packages. More in detail, the models coming from \code{rugarch} are: model Valid models (currently implemented) are 
#' 'sGARCH', 'eGARCH', 'gjrGARCH', 'iGARCH', and 'csGARCH'. The models implemented from \code{rumidas} are: 
#' 'GM_skew','GM_noskew', 'DAGM_skew', and 'DAGM_noskew'
#' @param distribution **optional** Distribution chosen for the univariate estimation. Valid choices are: "norm" (by default) and "std", 
#' respectively, for the Normal and Student's t distributions
#' @param MV **optional** MIDAS variable to include in the univariate estimation, if the model specificied is a GARCH-MIDAS 
#' (GM, \insertCite{engle_ghysels_sohn_2013;textual}{dccmidas}) or a Double Asymmetric GM (DAGM, \insertCite{engle_ghysels_sohn_2013;textual}{dccmidas}). 
#' In the case of MIDAS-based models, please provide a list of the MIDAS variables obtained from the \link[rumidas]{mv_into_mat} function. 
#' If the same MV variable is used,  then provide always a list, with the same (transformed) variable repeated
#' @param K **optional** The number of lagged realization of MV variable to use, if 'univ_model' has a MIDAS term
#' @param corr_model Correlation model used. Valid choices are: "cDCC" (the corrected DCC of \insertCite{aielli2013dynamic;textual}{dccmidas}),
#' "aDCC" (the asymmetric DCC model of \insertCite{cappiello2006asymmetric;textual}{dccmidas}),
#' "DECO" (Dynamic equicorrelation of \insertCite{engle2012dynamic;textual}{dccmidas}), 
#' and "DCCMIDAS" (the DCC-MIDAS of \insertCite{colacito2011component;textual}{dccmidas}). By detault, it is "cDCC" 
#' @param lag_fun **optional**. Lag function to use. Valid choices are "Beta" (by default) and "Almon", 
#' for the Beta and Exponential Almon lag functions, respectively, if 'univ_model' has a 
#' MIDAS term and/or if 'corr_model' is "DCCMIDAS"
#' @param N_c **optional** Number of (lagged) realizations to use for the standarized residuals forming the long-run correlation, if 'corr_model' is "DCCMIDAS"
#' @param K_c **optional** Number of (lagged) realizations to use for the long-run correlation, if 'corr_model' is "DCCMIDAS"
#' @param out_of_sample **optional** A positive integer indicating the number of periods before the last to keep for out of sample forecasting
#' @return \code{dcc_fit} returns an object of class 'dccmidas'. The function \code{\link{summary.dccmidas}} 
#' can be used to print a summary of the results. Moreover, an object of class 'dccmidas' is a list containing the following components:
#' \itemize{
#' 	\item assets: Names of the assets considered.
#'   \item model: Univariate model used in the first step.
#'    \item est_univ_model: List of matrixes of estimated coefficients of the univariate model, with the QML \insertCite{Bollerslev_Wooldridge_1992}{dccmidas} standard errors. 
#'    \item corr_coef_mat: Matrix of estimated coefficients of the correlation model, with the QML standard errors.
#'    \item mult_model: Correlation model used in the second step.
#'   \item obs: The number of daily observations used for the estimation.
#'   \item period: The period of the (in-sample) estimation.
#'   \item H_t: Conditional covariance matrix, reported as an array.
#'	\item R_t: Conditional correlation matrix, reported as an array.
#'	\item R_t_bar: Conditional long-run correlation matrix, reported as an array, if the correlation matrix includes a MIDAS specification.
#'   \item H_t_oos: Conditional covariance matrix, reported as an array, for the out-of-sample period, if present.
#'	\item R_t_oos: Conditional correlation matrix, reported as an array, for the out-of-sample period, if present.
#'	\item R_t_bar_oos: Conditional long-run correlation matrix, reported as an array, if the correlation matrix includes a MIDAS specification, for the out-of-sample period, if present.
#'   \item est_time: Time of estimation.
#'    \item Days: Days of the (in-)sample period.
#'    \item llk: The value of the log-likelihood (for the second step) at the maximum.
#' }
#' @importFrom Rdpack reprompt
#' @import roll
#' @import rumidas
#' @import rugarch
#' @details
#' Function \code{dcc_fit} implements the two-steps estimation of the DCC models. In the first step, a variety of univariate models are
#' considered. These models can be selected using for the parameter 'univ_model' one of the following choices: 'sGARCH' 
#' (standard GARCH of \insertCite{bollerslev1986generalized;textual}{dccmidas}), 
#' 'eGARCH' of \insertCite{nelson1991conditional;textual}{dccmidas}, 
#' 'gjrGARCH' of \insertCite{glosten_1993;textual}{dccmidas}, 
#' 'iGARCH' (Integrated GARCH of \insertCite{engle1986modelling;textual}{dccmidas}), 
#' 'csGARCH' (the Component GARCH of \insertCite{Engle_lee_1999;textual}{dccmidas}),
#' 'GM_noskew' and 'GM_skew' (the GARCH-MIDAS model of \insertCite{engle_ghysels_sohn_2013;textual}{dccmidas}, respectively, 
#' without and with the asymmetric term in the short-run component),
#'  and 'DAGM_noskew' and 'DAGM_skew' (the Double Asymmetric GARCH-MIDAS model of \insertCite{amendola_candila_gallo_2019;textual}{dccmidas},
#' respectively, without and with the asymmetric term in the short-run component).
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' require(xts)
#' # daily log-returns
#' # close to close daily log-returns
#' r_t_s<-diff(log(sp500['2010/2019'][,3]))
#' r_t_s[1]<-0
#' r_t_n<-diff(log(nasdaq['2010/2019'][,3]))
#' r_t_n[1]<-0
#' r_t_f<-diff(log(ftse100['2010/2019'][,3]))
#' r_t_f[1]<-0
#' db_m<-merge.xts(r_t_s,r_t_n,r_t_f)
#' db_m<-db_m[complete.cases(db_m),]
#' colnames(db_m)<-c("S&P500","NASDAQ","FTSE100")
#' # list of returns
#' r_t<-list(db_m[,1],db_m[,2],db_m[,3])
#' # MV transformation (same MV for all the stocks)
#' require(rumidas)
#' mv_m<-mv_into_mat(r_t[[1]],diff(indpro),K=12,"monthly")
#' # list of MV
#' MV<-list(mv_m,mv_m,mv_m)
#' # estimation
#' K_c<-144
#' N_c<-36
#' dccmidas_est<-dcc_fit(r_t,univ_model="GM_noskew",distribution="norm",
#' MV=MV,K=12,corr_model="DCCMIDAS",N_c=N_c,K_c=K_c)
#' dccmidas_est
#' summary.dccmidas(dccmidas_est)
#' }
#' @export

dcc_fit<-function(r_t,univ_model="sGARCH",distribution="norm",
MV=NULL,K=NULL,corr_model="cDCC",lag_fun="Beta",N_c=NULL,K_c=NULL,out_of_sample=NULL){

################################### checks

cond_r_t<- class(r_t)

if(cond_r_t != "list") { stop(
cat("#Warning:\n Parameter 'r_t' must be a list object. Please provide it in the correct form \n")
)}

if ((univ_model %in% c("GM_noskew", "GM_skew", "DAGM_noskew", "DAGM_skew")) & 
(!inherits(MV, "list") | length(MV) != length(r_t) | is.null(K)))
{ stop(
cat("#Warning:\n If you want to estimate a GARCH-MIDAS model, please provide parameters MV and K in the correct form \n")
)}

if(
(corr_model=="DCCMIDAS"|corr_model=="ADCCMIDAS")&(is.null(N_c)|is.null(K_c))
) { stop(
cat("#Warning:\n If you want to estimate a DCC-MIDAS model, please provide parameters N_c and K_c \n")
)}

################################### merge (eventually)

Num_assets<-length(r_t)

len<-rep(NA,Num_assets)

for(i in 1:Num_assets){
len[i]<-length(r_t[[i]])
}

if(stats::var(len)!=0){
db<- do.call(xts::merge.xts,r_t)
db<-db[stats::complete.cases(db),] } else
{
db<-do.call(xts::merge.xts,r_t)
}

for(i in 1:Num_assets){
colnames(db)[i]<-colnames(r_t[[i]])
}


################################### in-sample and out-of-sample

TT<-nrow(db)

if (missing(out_of_sample)){ # out_of_sample parameter is missing

sd_m<-matrix(NA,ncol=Num_assets,nrow=TT)

estim_period<-range(stats::time(db))
days_period<-stats::time(db)
est_obs<-nrow(db)

} else {					# out_of_sample parameter is present

sd_m<-matrix(NA,ncol=Num_assets,nrow=TT-out_of_sample)
db_oos<-db[(TT-out_of_sample+1):TT,]
sd_m_oos<-matrix(NA,ncol=Num_assets,nrow=out_of_sample)

estim_period<-range(stats::time(db[1:(TT-out_of_sample),]))
days_period<-stats::time(db[1:(TT-out_of_sample),])
est_obs<-nrow(db[1:(TT-out_of_sample),])

}


################################### setting the estimation list

est_details<-list()

likelihood_at_max<-list()

start<-Sys.time()

########################################### first step (univariate GARCH)

if(univ_model=="sGARCH"){

uspec<-rugarch::ugarchspec(
variance.model=list(model="sGARCH", garchOrder=c(1,1)),
mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
distribution.model = distribution)


if (missing(out_of_sample)){

for(i in 1:Num_assets){
u_est<-rugarch::ugarchfit(spec=uspec,data=db[,i])
est_details[[i]]<-u_est@fit$robust.matcoef
likelihood_at_max[[i]]<-u_est@fit$LLH
sd_m[,i]<-u_est@fit$sigma
}

} else {

for(i in 1:Num_assets){
u_est<-rugarch::ugarchfit(spec=uspec,data=db[,i],out.sample=out_of_sample)
est_details[[i]]<-u_est@fit$robust.matcoef
likelihood_at_max[[i]]<-u_est@fit$LLH
sd_m[,i]<-u_est@fit$sigma
sd_m_oos[,i]<-as.numeric(sigma(rugarch::ugarchforecast(u_est,  n.ahead = 1, 
n.roll = c(out_of_sample-1), 
out.sample = c(out_of_sample-1))))
}


} 
} else if (univ_model=="gjrGARCH"){

uspec<-rugarch::ugarchspec(
variance.model=list(model="gjrGARCH", garchOrder=c(1,1)),
mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
distribution.model = distribution)

if (missing(out_of_sample)){

for(i in 1:Num_assets){
u_est<-rugarch::ugarchfit(spec=uspec,data=db[,i])
est_details[[i]]<-u_est@fit$robust.matcoef
likelihood_at_max[[i]]<-u_est@fit$LLH
sd_m[,i]<-u_est@fit$sigma
}

} else {

for(i in 1:Num_assets){
u_est<-rugarch::ugarchfit(spec=uspec,data=db[,i],out.sample=out_of_sample)
est_details[[i]]<-u_est@fit$robust.matcoef
likelihood_at_max[[i]]<-u_est@fit$LLH
sd_m[,i]<-u_est@fit$sigma
sd_m_oos[,i]<-as.numeric(sigma(rugarch::ugarchforecast(u_est,  n.ahead = 1, 
n.roll = c(out_of_sample-1), 
out.sample = c(out_of_sample-1))))
}


}
} else if (univ_model=="eGARCH"){

uspec<-rugarch::ugarchspec(
variance.model=list(model="eGARCH", garchOrder=c(1,1)),
mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
distribution.model = distribution)

if (missing(out_of_sample)){

for(i in 1:Num_assets){
u_est<-rugarch::ugarchfit(spec=uspec,data=db[,i])
est_details[[i]]<-u_est@fit$robust.matcoef
likelihood_at_max[[i]]<-u_est@fit$LLH
sd_m[,i]<-u_est@fit$sigma
}

} else {

for(i in 1:Num_assets){
u_est<-rugarch::ugarchfit(spec=uspec,data=db[,i],out.sample=out_of_sample)
est_details[[i]]<-u_est@fit$robust.matcoef
likelihood_at_max[[i]]<-u_est@fit$LLH
sd_m[,i]<-u_est@fit$sigma
sd_m_oos[,i]<-as.numeric(sigma(rugarch::ugarchforecast(u_est,  n.ahead = 1, 
n.roll = c(out_of_sample-1), 
out.sample = c(out_of_sample-1))))
}


}
} else if (univ_model=="iGARCH"){

uspec<-rugarch::ugarchspec(
variance.model=list(model="iGARCH", garchOrder=c(1,1)),
mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
distribution.model = distribution)

if (missing(out_of_sample)){

for(i in 1:Num_assets){
u_est<-rugarch::ugarchfit(spec=uspec,data=db[,i])
est_details[[i]]<-u_est@fit$robust.matcoef
likelihood_at_max[[i]]<-u_est@fit$LLH
sd_m[,i]<-u_est@fit$sigma
}

} else {

for(i in 1:Num_assets){
u_est<-rugarch::ugarchfit(spec=uspec,data=db[,i],out.sample=out_of_sample)
est_details[[i]]<-u_est@fit$robust.matcoef
likelihood_at_max[[i]]<-u_est@fit$LLH
sd_m[,i]<-u_est@fit$sigma
sd_m_oos[,i]<-as.numeric(sigma(rugarch::ugarchforecast(u_est,  n.ahead = 1, 
n.roll = c(out_of_sample-1), 
out.sample = c(out_of_sample-1))))
}


}
} else if (univ_model=="csGARCH"){

uspec<-rugarch::ugarchspec(
variance.model=list(model="csGARCH", garchOrder=c(1,1)),
mean.model=list(armaOrder=c(0,0), include.mean=FALSE),
distribution.model = distribution)

if (missing(out_of_sample)){

for(i in 1:Num_assets){
u_est<-rugarch::ugarchfit(spec=uspec,data=db[,i])
est_details[[i]]<-u_est@fit$robust.matcoef
likelihood_at_max[[i]]<-u_est@fit$LLH
sd_m[,i]<-u_est@fit$sigma
}

} else {

for(i in 1:Num_assets){
u_est<-rugarch::ugarchfit(spec=uspec,data=db[,i],out.sample=out_of_sample)
est_details[[i]]<-u_est@fit$robust.matcoef
likelihood_at_max[[i]]<-u_est@fit$LLH
sd_m[,i]<-u_est@fit$sigma
sd_m_oos[,i]<-as.numeric(sigma(rugarch::ugarchforecast(u_est,  n.ahead = 1, 
n.roll = c(out_of_sample-1), 
out.sample = c(out_of_sample-1))))
}


}
} else if (univ_model=="GM_noskew"){ 

if (missing(out_of_sample)){

for(i in 1:Num_assets){
u_est<-rumidas::ugmfit(model="GM",skew="NO",distribution=distribution,db[,i],MV[[i]],K=K,lag_fun=lag_fun)
est_details[[i]]<-u_est$rob_coef_mat
likelihood_at_max[[i]]<-u_est$loglik
sd_m[,i]<-u_est$est_vol_in_s

} 

} else {

for(i in 1:Num_assets){
u_est<-rumidas::ugmfit(model="GM",skew="NO",distribution=distribution,db[,i],MV[[i]],K=K,lag_fun=lag_fun,out_of_sample=out_of_sample)
est_details[[i]]<-u_est$rob_coef_mat
likelihood_at_max[[i]]<-u_est$loglik
sd_m[,i]<-u_est$est_vol_in_s
sd_m_oos[,i]<-zoo::coredata(u_est$est_vol_oos)
} 
} 
} else if (univ_model=="GM_skew"){ 

if (missing(out_of_sample)){

for(i in 1:Num_assets){
u_est<-rumidas::ugmfit(model="GM",skew="YES",distribution=distribution,db[,i],MV[[i]],K=K,lag_fun=lag_fun)
est_details[[i]]<-u_est$rob_coef_mat
likelihood_at_max[[i]]<-u_est$loglik
sd_m[,i]<-u_est$est_vol_in_s

} 

} else {

for(i in 1:Num_assets){
u_est<-rumidas::ugmfit(model="GM",skew="NO",distribution=distribution,db[,i],MV[[i]],K=K,lag_fun=lag_fun,out_of_sample=out_of_sample)
est_details[[i]]<-u_est$rob_coef_mat
likelihood_at_max[[i]]<-u_est$loglik
sd_m[,i]<-u_est$est_vol_in_s
sd_m_oos[,i]<-zoo::coredata(u_est$est_vol_oos)
} 
} 
 
} else if (univ_model=="DAGM_noskew"){ 

if (missing(out_of_sample)){

for(i in 1:Num_assets){
u_est<-rumidas::ugmfit(model="DAGM",skew="NO",distribution=distribution,db[,i],MV[[i]],K=K,lag_fun=lag_fun)
est_details[[i]]<-u_est$rob_coef_mat
likelihood_at_max[[i]]<-u_est$loglik
sd_m[,i]<-u_est$est_vol_in_s

} 

} else {

for(i in 1:Num_assets){
u_est<-rumidas::ugmfit(model="GM",skew="NO",distribution=distribution,db[,i],MV[[i]],K=K,lag_fun=lag_fun,out_of_sample=out_of_sample)
est_details[[i]]<-u_est$rob_coef_mat
likelihood_at_max[[i]]<-u_est$loglik
sd_m[,i]<-u_est$est_vol_in_s
sd_m_oos[,i]<-zoo::coredata(u_est$est_vol_oos)
} 
} 
} else if (univ_model=="DAGM_skew"){ 
if (missing(out_of_sample)){

for(i in 1:Num_assets){
u_est<-rumidas::ugmfit(model="DAGM",skew="YES",distribution=distribution,db[,i],MV[[i]],K=K,lag_fun=lag_fun)
est_details[[i]]<-u_est$rob_coef_mat
likelihood_at_max[[i]]<-u_est$loglik
sd_m[,i]<-u_est$est_vol_in_s

} 

} else {

for(i in 1:Num_assets){
u_est<-rumidas::ugmfit(model="GM",skew="NO",distribution=distribution,db[,i],MV[[i]],K=K,lag_fun=lag_fun,out_of_sample=out_of_sample)
est_details[[i]]<-u_est$rob_coef_mat
likelihood_at_max[[i]]<-u_est$loglik
sd_m[,i]<-u_est$est_vol_in_s
sd_m_oos[,i]<-zoo::coredata(u_est$est_vol_oos)
} 
} 
}

cat("First step: completed \n")

##### standardized residuals

if (missing(out_of_sample)){

D_t<-array(0,dim=c(Num_assets,Num_assets,TT))
eps_t<-array(0,dim=c(Num_assets,1,TT))
db_a<-array(NA,dim=c(Num_assets,1,TT))
db_no_xts<-zoo::coredata(db)

for(tt in 1:TT){
db_a[,,tt]<-t(db_no_xts[tt,])
diag(D_t[,,tt])<-sd_m[tt,]
eps_t[,,tt]<-Inv(D_t[,,tt])%*%db_a[,,tt]
}
} else {

## in-sample
D_t<-array(0,dim=c(Num_assets,Num_assets,TT-out_of_sample))
eps_t<-array(0,dim=c(Num_assets,1,TT-out_of_sample))
db_a<-array(NA,dim=c(Num_assets,1,TT-out_of_sample))
db_no_xts<-zoo::coredata(db[1:(TT-out_of_sample),])

for(tt in 1:(TT-out_of_sample)){
db_a[,,tt]<-t(db_no_xts[tt,])
diag(D_t[,,tt])<-sd_m[tt,]
eps_t[,,tt]<-Inv(D_t[,,tt])%*%db_a[,,tt]
}

## out-of-sample
D_t_oos<-array(0,dim=c(Num_assets,Num_assets,out_of_sample))
eps_t_oos<-array(0,dim=c(Num_assets,1,out_of_sample))
db_a_oos<-array(NA,dim=c(Num_assets,1,out_of_sample))
db_no_xts_oos<-zoo::coredata(db_oos)

for(tt in 1:out_of_sample){
db_a_oos[,,tt]<-t(db_no_xts_oos[tt,])
diag(D_t_oos[,,tt])<-sd_m_oos[tt,]
eps_t_oos[,,tt]<-Inv(D_t_oos[,,tt])%*%db_a_oos[,,tt]
}


}


########################################### second step estimation (cDCC, aDCC, DECO or DCC-MIDAS)

if (corr_model=="DCCMIDAS"&lag_fun=="Beta"){

start_val<-ui<-ci<-NULL

start_val<-c(a=0.01,b=0.8,w2=2)

ui<-rbind(
c(1,0,0),       	 	 	## alpha>0.0001
c(0,1,0),        		 	## beta>0.001
c(-1,-1,0),	     		## alpha+beta<1
c(0,0,1))				## w2>1.001

ci<-c(-0.0001,-0.001,0.999,-1.001)

m_est<-maxLik(logLik=dccmidas_loglik,
start=start_val,
res=eps_t,
lag_fun=lag_fun,
N_c=N_c,
K_c=K_c,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS") 

} else if (corr_model=="DCCMIDAS"&lag_fun=="Almon") {

start_val<-ui<-ci<-NULL

start_val<-c(a=0.01,b=0.8,w2=-0.1)

ui<-rbind(
c(1,0,0),       	 	 	## alpha>0.0001
c(0,1,0),        		 	## beta>0.001
c(-1,-1,0),	     		## alpha+beta<1
c(0,0,-1))					## w2<0

ci<-c(-0.0001,-0.001,0.999,0)

m_est<-maxLik(logLik=dccmidas_loglik,
start=start_val,
res=eps_t,
lag_fun=lag_fun,
N_c=N_c,
K_c=K_c,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS") 

} else if (corr_model=="ADCCMIDAS"&lag_fun=="Beta"){

start_val<-ui<-ci<-NULL

start_val<-c(a=0.01,b=0.8,g=0.01,w2=2)

eta<-array(0,dim=c(Num_assets,1,TT))
eta<-ifelse(eps_t<0,1,0)*eps_t

sample_S<-stats::cov(t(apply(eps_t, 3L, c)))
sample_N<-stats::cov(t(apply(eta, 3L, c)))

inv_sample_S<-sample_S%^%(-0.5)
inv_sample_N<-Inv(sample_N)

delta<-eigen(inv_sample_S%*%inv_sample_N%*%inv_sample_S)$values[1]

ui<-rbind(
c(1,0,0,0),       	 	 ## a>0.0001
c(0,1,0,0),        		 ## b>0.001
c(0,0,1,0),				 ## g>0.001
c(-1,-1,-delta,0),	     ## a+b+delta*g<1
c(0,0,-1,0),				 ## g<0.05
c(0,0,0,1))				 ## w2>1.001

ci<-c(-0.0001,-0.001,-0.001,0.999,0.05,-1.001)

m_est<-maxLik(logLik=a_dccmidas_loglik,
start=start_val,
res=eps_t,
lag_fun=lag_fun,
N_c=N_c,
K_c=K_c,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS") 

} else if (corr_model=="ADCCMIDAS"&lag_fun=="Almon"){

start_val<-ui<-ci<-NULL

start_val<-c(a=0.01,b=0.8,g=0.01,w2=-0.1)


eta<-array(0,dim=c(Num_assets,1,TT))
eta<-ifelse(eps_t<0,1,0)*eps_t

sample_S<-stats::cov(t(apply(eps_t, 3L, c)))
sample_N<-stats::cov(t(apply(eta, 3L, c)))

inv_sample_S<-sample_S%^%(-0.5)
inv_sample_N<-Inv(sample_N)

delta<-eigen(inv_sample_S%*%inv_sample_N%*%inv_sample_S)$values[1]

ui<-rbind(
c(1,0,0,0),       	 	 ## a>0.0001
c(0,1,0,0),        		 ## b>0.001
c(0,0,1,0),				 ## g>0.001
c(-1,-1,-delta,0),	     ## a+b+delta*g<1
c(0,0,-1,0),				 ## g<0.05
c(0,0,0,-1))				 ## w2>1.001

ci<-c(-0.0001,-0.001,-0.001,0.999,0.05,0)

m_est<-maxLik(logLik=a_dccmidas_loglik,
start=start_val,
res=eps_t,
lag_fun=lag_fun,
N_c=N_c,
K_c=K_c,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS") 

} else if (corr_model=="cDCC"){

start_val<-ui<-ci<-NULL

start_val<-c(a=0.01,b=0.8)

ui<-rbind(
c(1,0),       	 	 	## alpha>0.0001
c(0,1),        		 	## beta>0.001
c(-1,-1))	     		## alpha+beta<1

ci<-c(-0.0001,-0.001,0.999)

m_est<-maxLik(logLik=dcc_loglik,
start=start_val,
res=eps_t,
K_c=K_c,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS") 

} else if (corr_model=="DECO"){

start_val<-ui<-ci<-NULL

start_val<-c(a=0.01,b=0.8)

ui<-rbind(
c(1,0),       	 	 	## alpha>0.0001
c(0,1),        		 	## beta>0.001
c(-1,-1))	     		## alpha+beta<1

ci<-c(-0.0001,-0.001,0.999)

m_est<-maxLik(logLik=deco_loglik,
start=start_val,
res=eps_t,
K_c=K_c,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS") 

} else if (corr_model=="aDCC"){

start_val<-ui<-ci<-NULL

start_val<-c(a=0.01,b=0.8,g=0.1)

eta<-array(0,dim=c(Num_assets,1,TT))
eta<-ifelse(eps_t<0,1,0)*eps_t

sample_S<-stats::cov(t(apply(eps_t, 3L, c)))
sample_N<-stats::cov(t(apply(eta, 3L, c)))

inv_sample_S<-sample_S%^%(-0.5)
inv_sample_N<-Inv(sample_N)

delta<-eigen(inv_sample_S%*%inv_sample_N%*%inv_sample_S)$values[1]

ui<-rbind(
c(1,0,0),       	 	 	## a>0.0001
c(0,1,0),        		 	## b>0.001
c(0,0,1),				## g>0.001
c(-1,-1,-delta),	     	## a+b+delta*g<1
c(0,0,-1))				## g<0.15

ci<-c(-0.0001,-0.001,-0.001,0.999,0.15)

m_est<-maxLik(logLik=a_dcc_loglik,
start=start_val,
res=eps_t,
K_c=K_c,
constraints=list(ineqA=ui, ineqB=ci),
iterlim=1000,
method="BFGS") 

}

########################################### end second step

###### matrix of coefficients (second step)

est_coef<-stats::coef(m_est)
N_coef<-length(est_coef)

mat_coef<-data.frame(rep(NA,N_coef),rep(NA,N_coef),rep(NA,N_coef),rep(NA,N_coef))
colnames(mat_coef)<-c("Estimate","Std. Error","t value","Pr(>|t|)")

rownames(mat_coef)<-names(est_coef)

mat_coef[,1]<-round(est_coef,6)
mat_coef[,2]<-round(QMLE_sd(m_est),6)
mat_coef[,3]<-round(est_coef/QMLE_sd(m_est),6)
mat_coef[,4]<-round(apply(rbind(est_coef/QMLE_sd(m_est)),1,function(x) 2*(1-stats::pnorm(abs(x)))),6)

if(corr_model=="DCCMIDAS"){

dcc_mat_est_fin<-dccmidas_mat_est(est_coef,eps_t,D_t,lag_fun=lag_fun,N_c=N_c,K_c=K_c)

if (!missing(out_of_sample)){

dcc_mat_est_fin_oos<-dccmidas_mat_est(est_coef,eps_t_oos,D_t_oos,
lag_fun=lag_fun,N_c=N_c,K_c=K_c)

} else {

dcc_mat_est_fin_oos<-list(NA,NA,NA)

}


} else if (corr_model=="cDCC") {

dcc_mat_est_fin<-dcc_mat_est(est_coef,eps_t,D_t,K_c=K_c)


if (!missing(out_of_sample)){

dcc_mat_est_fin_oos<-dcc_mat_est(est_coef,
eps_t_oos,D_t_oos,K_c=K_c)

} else {

dcc_mat_est_fin_oos<-list(NA,NA)

}


} else if (corr_model=="aDCC") {

dcc_mat_est_fin<-a_dcc_mat_est(est_coef,eps_t,D_t,K_c=K_c)


if (!missing(out_of_sample)){

dcc_mat_est_fin_oos<-a_dcc_mat_est(est_coef,
eps_t_oos,D_t_oos,K_c=K_c)

} else {

dcc_mat_est_fin_oos<-list(NA,NA)

}



} else if (corr_model=="ADCCMIDAS") {

dcc_mat_est_fin<-a_dccmidas_mat_est(est_coef,eps_t,D_t,lag_fun=lag_fun,N_c=N_c,K_c=K_c)

if (!missing(out_of_sample)){

dcc_mat_est_fin_oos<-a_dccmidas_mat_est(est_coef,
eps_t_oos,D_t_oos,lag_fun=lag_fun,N_c=N_c,K_c=K_c)

} else {

dcc_mat_est_fin_oos<-list(NA,NA,NA)

}


} else if (corr_model=="DECO") {

dcc_mat_est_fin<-deco_mat_est(est_coef,eps_t,D_t,K_c=K_c)


if (!missing(out_of_sample)){

dcc_mat_est_fin_oos<-deco_mat_est(est_coef,
eps_t_oos,D_t_oos,K_c=K_c)

} else {

dcc_mat_est_fin_oos<-list(NA,NA)

}

}

######

end<-Sys.time()-start


if(corr_model=="DCCMIDAS"|corr_model=="ADCCMIDAS"){
fin_res<-list(
assets=colnames(db),
model=univ_model,
est_univ_model=est_details,
corr_coef_mat=mat_coef,
mult_model=corr_model,
obs=est_obs,
period=estim_period,
"H_t"=dcc_mat_est_fin[[1]],
"R_t"=dcc_mat_est_fin[[2]],
"R_t_bar"=dcc_mat_est_fin[[3]],
"H_t_oos"=dcc_mat_est_fin_oos[[1]],
"R_t_oos"=dcc_mat_est_fin_oos[[2]],
"R_t_bar_oos"=dcc_mat_est_fin_oos[[3]],
est_time=end,
Days=days_period,
llk=stats::logLik(m_est)
)
} else {
fin_res<-list(
assets=colnames(db),
model=univ_model,
est_univ_model=est_details,
corr_coef_mat=mat_coef,
mult_model=corr_model,
obs=est_obs,
period=estim_period,
"H_t"=dcc_mat_est_fin[[1]],
"R_t"=dcc_mat_est_fin[[2]],
"H_t_oos"=dcc_mat_est_fin_oos[[1]],
"R_t_oos"=dcc_mat_est_fin_oos[[2]],
est_time=end,
Days=days_period,
llk=stats::logLik(m_est)
)
}

class(fin_res)<-c("dccmidas")
return(fin_res)
print.dccmidas(fin_res)


}

#' Var-cov matrix evaluation 
#'
#' Evaluates the estimated var-cov matrix H_t with respect to a covariance proxy, under different robust loss functions 
#' \insertCite{laurent2013loss}{dccmidas}. The losses considered are also used in \insertCite{amendola_2020;textual}{dccmidas}.
#' @param H_t Estimated covariance matrix, formatted as array
#' @param cov_proxy **optional** Covariance matrix, formatted as array
#' @param r_t **optional** List of daily returns used to calculate H_t. If parameter 'cov_proxy' is not provided, then
#' r_t must be included. In this case, a (noise) proxy will be automatically used
#' @param loss Robust loss function to use. Valid choices are: "FROB" for Frobenius (by default), "SFROB" for Squared Frobenius,
#' "EUCL" for Euclidean, "QLIKE" for QLIKE and "RMSE" for Root Mean Squared Errors
#' @return The value of the loss for each \eqn{t}
#' @importFrom Rdpack reprompt
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' require(xts)
#' # close to close daily log-returns
#' r_t_s<-diff(log(sp500['2010/2019'][,3]))
#' r_t_s[1]<-0
#' r_t_n<-diff(log(nasdaq['2010/2019'][,3]))
#' r_t_n[1]<-0
#' r_t_f<-diff(log(ftse100['2010/2019'][,3]))
#' r_t_f[1]<-0
#' db_m<-merge.xts(r_t_s,r_t_n,r_t_f)
#' db_m<-db_m[complete.cases(db_m),]
#' colnames(db_m)<-c("S&P500","NASDAQ","FTSE100")
#' # list of returns
#' r_t<-list(db_m[,1],db_m[,2],db_m[,3])
#' # estimation
#' K_c<-144
#' N_c<-36
#' cdcc_est<-dcc_fit(r_t,univ_model="sGARCH",distribution="norm",
#' corr_model="DCCMIDAS",N_c=N_c,K_c=K_c)
#' cov_eval(cdcc_est$H_t,r_t=r_t)[(K_c+1):dim(cdcc_est$H_t)[3]]
#' }
#' @export

cov_eval<-function(H_t,cov_proxy=NULL,r_t=NULL,loss="FROB"){

################################### checks

if (is.null(cov_proxy)&is.null(r_t)){ stop(
cat("#Warning:\n At least one parameter between 'cov_proxy' and 'r_t' must be provided \n")
)}

cond_r_t<- class(r_t)

if(!is.null(r_t)&cond_r_t != "list") { stop(
cat("#Warning:\n Parameter 'r_t' must be a list object. Please provide it in the correct form \n")
)}


if(!is.null(cov_proxy)&!inherits(cov_proxy, "array")) { stop(
cat("#Warning:\n Parameter 'cov_proxy' must be an array object. Please provide it in the correct form \n")
)}


################################### merge (eventually)

if (!is.null(r_t)){

Num_assets<-length(r_t)

len<-rep(NA,Num_assets)

for(i in 1:Num_assets){
len[i]<-length(r_t[[i]])
}

if(stats::var(len)!=0){
db<- do.call(xts::merge.xts,r_t)
db<-db[stats::complete.cases(db),] 
} else {
db<-do.call(xts::merge.xts,r_t)
}

for(i in 1:Num_assets){
colnames(db)[i]<-colnames(r_t[[i]])
}

TT<-nrow(db)

} else {

Num_assets<-dim(cov_proxy)[2]

TT<-dim(cov_proxy)[3]

}

################################### calculate the cov matrix from the r_t array

if (!is.null(r_t)){

cov_proxy<-array(0,dim=c(Num_assets,Num_assets,TT))

db<-zoo::coredata(db)

for (tt in 1:TT){
cov_proxy[,,tt]<-cbind(db[tt,])%*%rbind(db[tt,])
}

}

################################## perform the evaluation

loss_v<-rep(NA,TT)

if (loss=="FROB"){

for (tt in 1:TT){

dif<-(H_t[,,tt]-cov_proxy[,,tt])

loss_v[tt]<-sum(diag( t(dif)%*%dif ) )

}
} else if (loss=="SFROB"){

for (tt in 1:TT){

dif_sq<- (H_t[,,tt]-cov_proxy[,,tt])^2

loss_v[tt]<-sum(eigen(dif_sq)$values)

}
} else if (loss=="EUCL"){

for (tt in 1:TT){

h_t_vec<- rbind(H_t[,,tt][lower.tri(H_t[,,tt],diag=TRUE)])
cov_proxy_vec<-rbind(cov_proxy[,,tt][lower.tri(cov_proxy[,,tt],diag=TRUE)])

dif_vec<-cov_proxy_vec-h_t_vec

loss_v[tt]<-dif_vec%*%t(dif_vec)

}

} else if (loss=="QLIKE"){

for (tt in 1:TT){

log_est_det<-log(Det(H_t[,,tt]))
h_t_cov<-Inv(H_t[,,tt])%*%cov_proxy[,,tt]
h_t_vec<- rbind(h_t_cov[lower.tri(h_t_cov)])
h_t_vec<-h_t_vec%*%cbind(rep(1,length(h_t_vec)))

loss_v[tt]<-log_est_det+h_t_vec

} 

} else if (loss=="RMSE"){

for (tt in 1:TT){

dif<-(H_t[,,tt]-cov_proxy[,,tt])
dif_norm<-(norm(dif,type="F"))^0.5
loss_v[tt]<-dif_norm

} 

}

return(loss_v)

}

#' BEKK fit 
#'
#' Obtains the estimation the scalar and diagonal BEKK model
#' @param r_t List of daily returns. At the moment, at most 5 assets can be considered
#' @param model Valid choices are: 'sBEKK'(scalar BEKK) and 'dBEKK' (diagonal BEKK)
#' @param R **optional** Number of random samples drawn from a Uniform distribution used
#' to inizialize the log-likelihood. Equal to 100 by default
#' @param out_of_sample **optional** A positive integer indicating the number of periods before the last to keep for out of sample forecasting
#' @return \code{bekk_fit} returns a list containing the following components:
#' \itemize{
#' 	\item assets: Names of the assets considered.
#'   \item mat_coef: Matrix of estimated coefficients of the model, with the QML standard errors.
#'   \item obs: The number of daily observations used for the estimation.
#'   \item period: The period of the estimation.
#'   \item H_t: Conditional covariance matrix, reported as an array. It refers to the in-sample
#'   period.
#'   \item est_time: Time of estimation.
#'   \item llk: The value of the log-likelihood at the maximum.
#'   \item H_t_oos: Conditional covariance matrix, reported as an array, for the out-of-sample period,
#'  if the param 'out_of_sample' is used.
#'  \item Days: Days of the (in-)sample period.
#' }
#' @importFrom Rdpack reprompt
#' @importFrom stats runif
#' @details
#' Function \code{bekk_fit} implements the estimation of scalar and diagonal BEKK models. For details on BEKK models, see \insertCite{engle1995multivariate;textual}{dccmidas}
#' @references
#' \insertAllCited{} 
#' @examples
#' \donttest{
#' require(xts)
#' # close to close daily log-returns
#' r_t_s<-diff(log(sp500['2010/2019'][,3]))
#' r_t_s[1]<-0
#' r_t_n<-diff(log(nasdaq['2010/2019'][,3]))
#' r_t_n[1]<-0
#' r_t_f<-diff(log(ftse100['2010/2019'][,3]))
#' r_t_f[1]<-0
#' db_m<-merge.xts(r_t_s,r_t_n,r_t_f)
#' db_m<-db_m[complete.cases(db_m),]
#' colnames(db_m)<-c("S&P500","NASDAQ","FTSE100")
#' # list of returns
#' r_t<-list(db_m[,1],db_m[,2],db_m[,3])
#' bekk_est<-bekk_fit(r_t,model="sBEKK")
#' bekk_est$mat_coef
#' }
#' @export

bekk_fit<-function(r_t,model="sBEKK",R=100, out_of_sample = NULL){

################################### checks

cond_r_t<- class(r_t)

if(cond_r_t != "list") { stop(
cat("#Warning:\n Parameter 'r_t' must be a list object. Please provide it in the correct form \n")
)}

if(length(r_t) > 5) { stop(
cat("#Warning:\n Parameter 'r_t' has to have at most five assets.  \n")
)}

################################### merge (eventually)

Num_assets<-length(r_t)

len<-rep(NA,Num_assets)

for(i in 1:Num_assets){
len[i]<-length(r_t[[i]])
}


db<- do.call(xts::merge.xts,r_t)
db<-db[stats::complete.cases(db),] 

for(i in 1:Num_assets){
colnames(db)[i]<-colnames(r_t[[i]])
}

################################### in-sample and out-of-sample

TT<-nrow(db)

if (missing(out_of_sample)){ # out_of_sample parameter is missing

db_est<-db
estim_period<-range(stats::time(db_est))
days_period<-stats::time(db_est)

} else {					# out_of_sample parameter is present

db_est<-db[1:(TT-out_of_sample),]
db_oos<-db[(TT-out_of_sample+1):TT,]

estim_period<-range(stats::time(db_est))
days_period<-stats::time(db_est[1:(TT-out_of_sample),])


}

db<-db_est

################################### setting the estimation list

start<-Sys.time()

########################################### scalar BEKK

if(model=="sBEKK"){

start_val<-NULL

num_param<-ifelse(ncol(db)==2,5,ifelse(ncol(db)==3,8,ifelse(ncol(db)==4,12,17)))

start_val<-rep(0.1,num_param)

if (num_param==5){

names(start_val)<-c("c_11","c_21","c_22","a","b")

} else if (num_param==8){

names(start_val)<-c("c_11","c_21","c_22","c_31","c_32","c_33","a","b")

} else if (num_param==12) {

names(start_val)<-c("c_11","c_21","c_22","c_31","c_32","c_33",
"c_41","c_42","c_43","c_44","a","b")

} else if (num_param==17) {

names(start_val)<-c("c_11","c_21","c_22","c_31","c_32","c_33",
"c_41","c_42","c_43","c_44","c_51","c_52","c_53","c_54","c_55","a","b")

}

#### best starting values

mat_start_val<-matrix(stats::runif(num_param*R),ncol=R,nrow=num_param)
ll_start<-rep(NA,R)

for(i in 1:R){
ll_start[i]<-sum(sBEKK_loglik(mat_start_val[,i],zoo::coredata(db)))
}

start_val_f<-mat_start_val[,which.max(ll_start)]
names(start_val_f)<-names(start_val)

m_est<-maxLik(logLik=sBEKK_loglik,
start=start_val_f,
ret=zoo::coredata(db),
iterlim=1000,
method="BFGS") 

H_t_est<-sBEKK_mat_est(stats::coef(m_est),zoo::coredata(db))


if (missing(out_of_sample)){ # out_of_sample parameter is missing

H_t_oos<-NA

} else {					# out_of_sample parameter is present

H_t_oos<-sBEKK_mat_est(stats::coef(m_est),zoo::coredata(db_oos))

}



} else if (model=="dBEKK"){

start_val<-NULL

num_param<-ifelse(ncol(db)==2,7,ifelse(ncol(db)==3,12,ifelse(ncol(db)==4,18,25)))

start_val<-rep(0.5,num_param)

if (num_param==7){

names(start_val)<-c("c_11","c_21","c_22","a_11","a_22","b_11","b_22")

} else if (num_param==12){

names(start_val)<-c("c_11","c_21","c_22","c_31","c_32","c_33",
"a_11","a_22","a_33","b_11","b_22","b_33")

} else if (num_param==18) {

names(start_val)<-c("c_11","c_21","c_22","c_31","c_32","c_33",
"c_41","c_42","c_43","c_44",
"a_11","a_22","a_33","a_44","b_11","b_22","b_33","b_44")

} else if (num_param==25) {

names(start_val)<-c("c_11","c_21","c_22","c_31","c_32","c_33",
"c_41","c_42","c_43","c_44",
"c_51","c_52","c_53","c_54","c_55",
"a_11","a_22","a_33","a_44","a_55","b_11","b_22","b_33","b_44","b_55")

}

ui<-matrix(0,ncol=num_param,nrow=4)

ui[1,which(names(start_val)=="a_11")]<-1
ui[2,which(names(start_val)=="b_11")]<-1
ui[3,which(names(start_val)=="a_22")]<-1
ui[4,which(names(start_val)=="b_22")]<-1

ci<-c(-0.0001,-0.001,-0.0001,-0.001)

#### best starting values

mat_start_val<-matrix(stats::runif(num_param*R),ncol=R,nrow=num_param)
ll_start<-rep(NA,R)

for(i in 1:R){
ll_start[i]<-sum(dBEKK_loglik(mat_start_val[,i],zoo::coredata(db)))
}

start_val_f<-mat_start_val[,which.max(ll_start)]
names(start_val_f)<-names(start_val)


m_est<-suppressWarnings(maxLik(logLik=dBEKK_loglik,
start=start_val_f,
ret=zoo::coredata(db),
iterlim=1000,
constraints=list(ineqA=ui, ineqB=ci),
method="BFGS"))

H_t_est<-dBEKK_mat_est(stats::coef(m_est),zoo::coredata(db))

if (missing(out_of_sample)){ # out_of_sample parameter is missing

H_t_oos<-NA

} else {					# out_of_sample parameter is present

H_t_oos<-sBEKK_mat_est(stats::coef(m_est),zoo::coredata(db_oos))

}


} 

###### matrix of coefficients 

est_coef<-stats::coef(m_est)
N_coef<-length(est_coef)

mat_coef<-data.frame(rep(NA,N_coef),rep(NA,N_coef),rep(NA,N_coef),rep(NA,N_coef))
colnames(mat_coef)<-c("Estimate","Std. Error","t value","Pr(>|t|)")

rownames(mat_coef)<-names(est_coef)

mat_coef[,1]<-round(est_coef,6)
mat_coef[,2]<-round(QMLE_sd(m_est),6)
mat_coef[,3]<-round(est_coef/QMLE_sd(m_est),6)
mat_coef[,4]<-round(apply(rbind(est_coef/QMLE_sd(m_est)),1,function(x) 2*(1-stats::pnorm(abs(x)))),6)

######

end<-Sys.time()-start

fin_res<-list(
assets=colnames(db),
mat_coef=mat_coef,
obs=nrow(db),
period=estim_period,
H_t=H_t_est,
est_time=end,
llk=stats::logLik(m_est),
H_t_oos=H_t_oos,
Days=days_period
)

return(fin_res)

}


#' Plot method for 'dccmidas' class
#'
#' Plots of the conditional volatilities on the main diagonal and of the conditional correlations on the
#' extra-diagonal elements.
#' @param x An object of class 'dccmidas', that is the result of a call to \code{\link{dcc_fit}}.
#' @param K_c **optional** Number of (lagged) realizations to use for the long-run correlation, , if 'corr_model' is "DCCMIDAS"
#' @param vol_col **optional** Color of the volatility and correlation plots. "black" by default
#' @param long_run_col **optional** Color of the long-run correlation plots, if present. "red" by default
#' @param cex_axis **optional** Size of the x-axis. Default to 0.75
#' @param LWD **optional** Width of the plotted lines. Default to 2
#' @param asset_sub **optional** Numeric vector of selected assets to consider for the plot. NULL by default
#' @return No return value, called for side effects 
#' @export 
#' @examples
#' \donttest{
#' require(xts)
#' # close to close daily log-returns
#' r_t_s<-diff(log(sp500['2010/2019'][,3]))
#' r_t_s[1]<-0
#' r_t_n<-diff(log(nasdaq['2010/2019'][,3]))
#' r_t_n[1]<-0
#' r_t_f<-diff(log(ftse100['2010/2019'][,3]))
#' r_t_f[1]<-0
#' db_m<-merge.xts(r_t_s,r_t_n,r_t_f)
#' db_m<-db_m[complete.cases(db_m),]
#' colnames(db_m)<-c("S&P500","NASDAQ","FTSE100")
#' # list of returns
#' r_t<-list(db_m[,1],db_m[,2],db_m[,3])
#' # estimation
#' K_c<-144
#' N_c<-36
#' cdcc_est<-dcc_fit(r_t,univ_model="sGARCH",distribution="norm",
#' corr_model="DCCMIDAS",N_c=N_c,K_c=K_c)
#' plot_dccmidas(cdcc_est,K_c=144)
#' }

plot_dccmidas <- function(x,K_c=NULL,vol_col="black",long_run_col="red",
cex_axis=0.75,LWD=2,asset_sub=NULL){

################################### checks

class_x<-class(x)

if(class_x != "dccmidas") { stop(
cat("#Warning:\n Parameter 'x' must be a 'dccmidas' object. Please provide it in the correct form \n")
)}

if(is.null(K_c)){
K_c<-2
} else {
K_c<-K_c+1
}


################################### store objects

TT<-dim(x$H_t)[3]

if(is.null(asset_sub)){

H_t<-x$H_t[,,K_c:TT] #var-cov array
R_t<-x$R_t[,,K_c:TT] #corr array

# is long-run correlation present?

if(!is.null(x$R_t_bar)){
R_t_bar<-x$R_t_bar[,,K_c:TT]
} else {
R_t_bar<-H_t*0
}
} else {

H_t<-x$H_t[asset_sub,asset_sub,K_c:TT] #var-cov array
R_t<-x$R_t[asset_sub,asset_sub,K_c:TT] #corr array

# is long-run correlation present?

if(!is.null(x$R_t_bar)){
R_t_bar<-x$R_t_bar[asset_sub,asset_sub,K_c:TT]
} else {
R_t_bar<-H_t*0
}

}

# unify H_t and R_t

mat<-R_t 
corr_mat<-R_t_bar

for(tt in 1:dim(mat)[3]){
diag(mat[,,tt])<-diag(H_t[,,tt])^0.5
diag(corr_mat[,,tt])<-0
}

if(is.null(asset_sub)){
Num_col<-length(x$assets)
} else {
Num_col<-length(x$assets[asset_sub])
}

if(Num_col>=12&is.null(asset_sub)){ stop(
cat("#Warning:\n There are too many assets to build a readable plot. Please use the parameter 'asset_sub' \n")
)}

Days<-x$Days[K_c:TT]

#################################### time periods

if (length(Days)<250){

end_x<-xts::endpoints(Days,on="months")
end_x<-end_x+c(rep(1,length(end_x)-1),0)

tick_x<-c(seq(zoo::as.Date(Days[1]), zoo::as.Date(Days[length(Days)]), "months"),
zoo::as.Date(Days[length(Days)]))

tick_x<-format(tick_x,"%m/%Y")

} else {

end_x<-xts::endpoints(Days,on="years")
end_x<-end_x+c(rep(1,length(end_x)-1),0)

tick_x<-c(seq(zoo::as.Date(Days[1]), zoo::as.Date(Days[length(Days)]), "years"),
zoo::as.Date(Days[length(Days)]))

tick_x<-format(tick_x,"%m/%Y")

}
#################################### matrix id

matrix_id<-matrix(1:Num_col^2,ncol=Num_col,byrow=TRUE)
matrix_id_2<-data.frame(t(utils::combn(1:Num_col, 2)))
colnames(matrix_id_2)<-c("row","col")
matrix_id_2$Name<-NA

for(i in 1:nrow(matrix_id_2)){
matrix_id_2$Name[i]<-paste(
x$assets[[matrix_id_2$row[[i]]]],"-",
x$assets[[matrix_id_2$col[[i]]]]
)
}

if(is.null(asset_sub)){
assets_name<-sapply(x$assets, function (x) rep(x,Num_col))
} else {
assets_name<-sapply(x$assets[asset_sub], function (x) rep(x,Num_col))
}

matrix_id_2$id<-paste(matrix_id_2$row,"_",matrix_id_2$col,sep="")


######################################### plot

oldpar <- graphics::par(no.readonly = TRUE)    
on.exit(graphics::par(oldpar))            

graphics::par(mfrow = c(Num_col, Num_col),     
mai = c(0.3, 0.3, 0.3, 0.3))

case_to_cons<-matrix_id[upper.tri(matrix_id,diag=T)]

for(i in 1:Num_col^2){

# exclude cases
if(!i %in% case_to_cons){
graphics::plot.new()
} else { #including cases

row_chosen<-which(matrix_id==i, arr.ind=TRUE)[1]
col_chosen<-which(matrix_id==i, arr.ind=TRUE)[2]
row_col_chosen<-paste(row_chosen,"_",col_chosen,sep="")

mat_plot<-mat[row_chosen,col_chosen,]
corr_plot<-corr_mat[row_chosen,col_chosen,]

main_title<-ifelse(row_chosen==col_chosen,assets_name[i],
matrix_id_2$Name[matrix_id_2$id==row_col_chosen])

if(!all(corr_plot==0)){ #long-run correlation is present
y_lim_min<-min(mat_plot[1:length(mat_plot)],corr_plot[1:length(corr_plot)])
y_lim_max<-max(mat_plot[1:length(mat_plot)],corr_plot[1:length(corr_plot)])
y_lim<-c(y_lim_min,y_lim_max)
} else {
y_lim_min<-min(mat_plot[1:length(mat_plot)])
y_lim_max<-max(mat_plot[1:length(mat_plot)])
y_lim<-c(y_lim_min,y_lim_max)
}


plot(1:length(mat_plot),mat_plot,type="l",xaxt="n",
xlab="",ylab="",main = main_title, col=vol_col,lwd=LWD,ylim=y_lim,cex.axis=cex_axis)

if(!all(corr_plot==0)){
graphics::lines(1:length(corr_plot),corr_plot,col=long_run_col,
lwd=LWD,ylim=y_lim)
}

graphics::axis(side=1, 
end_x, 
tick_x,
cex.axis=cex_axis) 

graphics::grid(NA,NULL)
 
graphics::abline(v=end_x,h=NA,col="gray",lty=3)

}
}

}

