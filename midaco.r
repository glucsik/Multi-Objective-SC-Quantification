########################### GATEWAY HEADER #############################
#                           
#     _|      _|  _|_|_|  _|_|_|      _|_|      _|_|_|    _|_|    
#     _|_|  _|_|    _|    _|    _|  _|    _|  _|        _|    _|  
#     _|  _|  _|    _|    _|    _|  _|_|_|_|  _|        _|    _|  
#     _|      _|    _|    _|    _|  _|    _|  _|        _|    _|  
#     _|      _|  _|_|_|  _|_|_|    _|    _|    _|_|_|    _|_|  
#
#                                                   Version 6.0
#
########################################################################
#
#           See the MIDACO user manual for detailed information
#
########################################################################
#
#    Author (C) :   Dr. Martin Schlueter
#                   Information Initiative Center,
#                   Division of Large Scale Computing Systems,
#                   Hokkaido University, JAPAN.
#
#    Email :        info@midaco-solver.com
#
#    URL :          http://www.midaco-solver.com
#       
########################################################################
pacman::p_load(caret, dplyr, IRon, scam, parallel, doParallel, doMC, grwat, adc, hydroGOF, hydrostats, hydroEvents, tidyr, spatialEco, zoo, npreg, fANCOVA, forecast, FlowScreen, stinepack, imputeTS, missForest, Metrics, dplyr, tstools, tsbox, highfrequency, tidyverse, install = T)

## Load
df <- as.data.frame(read.csv("Sespe_Q_SC.csv"))

## Preallocate
df$Date <- as.Date(df$Date, "%m/%d/%Y")
df$time = 1:nrow(df)
strflow = as.numeric(df[,2])
baseq = matrix(strflow,ncol=5,length(strflow))
surfq = matrix(strflow,ncol=5,length(strflow))
wu = df[1:365,]
cd = df[(nrow(df)-364):nrow(df),]

## Impute and filter SC
# trim
df = tail(df,-365)
df = head(df,-365)

# Imputing missing SC values
imp  <- missForest(df[,-1], verbose = TRUE )
df[,2:4] <- imp$ximp

# Filter SC with LOESS
#SC_loess <- loess.as(df$time, df$SC, degree = 1, criterion = c("aicc", "gcv")[2], user.span = NULL, plot = T)
#df$SC = SC_loess$fitted

# add wu and cd for RDF
df = rbind(wu,df,cd)

midaco_run <- function( problem, option, key)
{
  
########################################################################
########################################################################
########################################################################

# cat('starting MIDACO\n')

o = problem.o 
n = problem.n 
ni= problem.ni
m = problem.m 
me= problem.me

xl = problem.xl
xu = problem.xu
x  = problem.x
  
maxeval = option.maxeval
maxtime = option.maxtime

printeval = option.printeval
save2file = option.save2file

param = double(length = 13)   

param[ 1] = option.param1 
param[ 2] = option.param2 
param[ 3] = option.param3 
param[ 4] = option.param4 
param[ 5] = option.param5 
param[ 6] = option.param6 
param[ 7] = option.param7 
param[ 8] = option.param8 
param[ 9] = option.param9 
param[10] = option.param10
param[11] = option.param11
param[12] = option.param12
param[13] = option.param13

P = option.parallel

if( P < 1 ) { P = 1 }

########################################################################
# Select printing mode: native R print commands or DLL based printing.
# Note: DLL based printing is slightly faster than native R printing, 
# but DLL based printing might not work on Windows platforms.
########################################################################
native_R_print = 1 # [ 1=YES / 0=NO ]
########################################################################
########################################################################
########################################################################  

if(.Platform$OS.type == "unix"){ 
      dyn.load("midacoR.so")  
} else { 
      dyn.load("midacoR.dll") 
}

########################################################################
midaco <- function(P,o,n,ni,m,me,x,f,g,xl,xu,iflag,istop,param,rw,lrw,iw,liw,pf,lpf,key) 
########################################################################
{ out <- .Fortran( "midaco", i01=as.integer(P),i02=as.integer(o),i03=as.integer(n),
      i04=as.integer(ni),i05=as.integer(m),i06=as.integer(me),i07=as.double(x),
      i08=as.double(f),i09=as.double(g),i10=as.double(xl),i11=as.double(xu),
      i12=as.integer(iflag),i13=as.integer(istop),i14=as.double(param),
      i15=as.double(rw),i16=as.integer(lrw),i17=as.integer(iw),
      i18 = as.integer(liw),i19=as.double(pf),i20=as.integer(lpf),i21 = as.character(key) )      
out <- list( out$i07, out$i08, out$i09, out$i12, out$i13, out$i15, out$i17, out$i19); return(out) }
########################################################################
midaco_print <- function(c,printeval,save2file,iflag,istop,f,g,x,xl,xu,
                         o,n,ni,m,me,rw,pf,maxeval,maxtime,param,P,T,key) 
########################################################################
{ out <- .Fortran( "midaco_print", i01=as.integer(c),i02=as.integer(printeval), 
      i03=as.integer(save2file),i04=as.integer(iflag),i05=as.integer(istop),  
      i06=as.double(f),i07=as.double(g),i08=as.double(x),i09=as.double(xl),
      i10=as.double(xu),i11=as.integer(o),i12=as.integer(n),i13=as.integer(ni),
      i14=as.integer(m),i15=as.integer(me),i16=as.double(rw),i17=as.double(pf),
      i18=as.integer(maxeval),i19=as.integer(maxtime),i20=as.double(param),
      i21=as.integer(P),i22=as.integer(T),i23=as.character(key) ); return( out$i04 ) }
########################################################################
########################################################################
########################################################################

fileA <- "MIDACO_SCREEN.TXT"
fileB <- "MIDACO_SOLUTION.TXT"
fileC <- "MIDACO_PARETOFRONT.TXT"
if (file.exists(fileA)) file.remove(fileA)
if (file.exists(fileB)) file.remove(fileB)
if (file.exists(fileC)) file.remove(fileC)

########################################################################
########################################################################
########################################################################
#     Initializations and Workspace Applocation
      f = double(length = o)
      g = double(length = m)
      lrw = 120*n+20*m+20*o+20*P+P*(m+2*o)+o*o+5000
      liw = 3*n+P+1000
      if ( o <= 1 ) { 
        lpf = 1
      } else {
        lpf = 1000 * (o+m+n) + 1 +100
        if ( param[10] >= 1) { lpf =  round(param[10]) * (o+m+n) + 1 }
        if ( param[10] <=-1) { lpf = -round(param[10]) * (o+m+n) + 1 }        
      }
      rw = double(length = lrw)
      iw = integer(length = liw)
      pf = double(length = lpf)
      iflag=0
      istop=0
      T = 0





if ( P <= 1 ) {
########################################################################
#         
#     Call MIDACO by Reverse Communication
#      
########################################################################
if( native_R_print == 0 ){
midaco_print(1,printeval,save2file,iflag,istop,f,g,x,xl,xu,
             o,n,ni,m,me,rw,pf,maxeval,maxtime,param,P,T,key)
}else{
midacoRprint(1,printeval,save2file,iflag,istop,f,g,x,xl,xu,
             o,n,ni,m,me,rw,pf,maxeval,maxtime,param,P,T,key)
}

while( istop == 0) 
{

      result = problem_function( f , g , x )

      f = result[[1]]
      g = result[[2]] 

      out = midaco(P,o,n,ni,m,me,x,f,g,xl,xu,iflag,istop,
                   param,rw,lrw,iw,liw,pf,lpf,key) 

      x     = out[[1]]
      f     = out[[2]]
      g     = out[[3]]
      iflag = out[[4]]
      istop = out[[5]]
      rw    = out[[6]]
      iw    = out[[7]]
      pf    = out[[8]]     

      if( is.na(iw[10]) )
      {
        # cat('catch and repair NA')
        iw[10] = 1000000
      } 

      if( native_R_print == 0 ){
      iflag = midaco_print(2,printeval,save2file,iflag,istop,f,g,x,xl,xu,
                           o,n,ni,m,me,rw,pf,maxeval,maxtime,param,P,T,key)
      }else{
      iflag = midacoRprint(2,printeval,save2file,iflag,istop,f,g,x,xl,xu,
                           o,n,ni,m,me,rw,pf,maxeval,maxtime,param,P,T,key)      
      }    
}
########################################################################
########################################################################
########################################################################
}






if ( P > 1 ) {
########################################################################
########################################################################
########################################################################

library(foreach)
#library(doMC)
#registerDoMC(P)
library(doSNOW)
cl<-makeCluster(P)
registerDoSNOW(cl)

fff = double(length = o*P)    
ggg = double(length = m*P) 
xxx = double(length = n*P) 
z = double(length = n) 
q = double(length = o) 
r = double(length = m) 

for(c in 1:P) 
{
    for(i in 1:n) 
    {
        xxx[(c-1)*n+i] = x[i] # starting point
    }
}
########################################################################
#         
#     Call MIDACO by Reverse Communication
#      
########################################################################
if( native_R_print == 0 ){
midaco_print(1,printeval,save2file,iflag,istop,f,g,x,xl,xu,
             o,n,ni,m,me,rw,pf,maxeval,maxtime,param,P,T,key)
}else{
midacoRprint(1,printeval,save2file,iflag,istop,f,g,x,xl,xu,
             o,n,ni,m,me,rw,pf,maxeval,maxtime,param,P,T,key)
}

while( istop == 0) 
{

      # loopoutput <- foreach(c=1:P,  .combine='c') %dopar% 
      loopoutput <-foreach(c=1:P,  .combine='c', .export = 'problem_function') %dopar%       
      {
            for(i in 1:n){ z[i] = xxx[(c-1)*n+i] }

            result = problem_function( q , r , z )

            q = result[[1]]
            r = result[[2]]

            loopoutput <- list( q , r )
      }     


      for(c in 1:P) 
      {
        r = loopoutput[[(c-1)*2+1]]
        q = loopoutput[[(c-1)*2+2]]

        for(i in 1:o){ fff[(c-1)*o+i] = r[i] }
        for(i in 1:m){ ggg[(c-1)*m+i] = q[i] }
      }

      out = midaco(P,o,n,ni,m,me,xxx,fff,ggg,xl,xu,iflag,istop,
                   param,rw,lrw,iw,liw,pf,lpf,key) 

      xxx   = out[[1]]
      fff   = out[[2]]
      ggg   = out[[3]]
      iflag = out[[4]]
      istop = out[[5]]
      rw    = out[[6]]
      iw    = out[[7]]
      pf    = out[[8]]       

      if( is.na(iw[10]) )
      {
        # cat('catch and repair NA')
        iw[10] = 1000000
      } 
      
      if( native_R_print == 0 ){
      iflag = midaco_print(2,printeval,save2file,iflag,istop,fff,ggg,xxx,xl,xu,
                           o,n,ni,m,me,rw,pf,maxeval,maxtime,param,P,T,key)
      }else{
      iflag = midacoRprint(2,printeval,save2file,iflag,istop,fff,ggg,xxx,xl,xu,
                           o,n,ni,m,me,rw,pf,maxeval,maxtime,param,P,T,key)
      }      
}

stopCluster(cl)

########################################################################
########################################################################
########################################################################
# Extract f,g,x from arrays
for(i in 1:o) { f[i] = fff[i] }
for(i in 1:m) { g[i] = ggg[i] }
for(i in 1:n) { x[i] = xxx[i] }
    
}
########################################################################
##########################  Output Arguments ###########################
########################################################################   

# cat('MIDACO finished\n')

output <- list( f, g, x, iflag )
return(output)
}
########################################################################
####################### END OF MIDACO R GATEWAY ########################
########################################################################














########################################################################
########################################################################
########################################################################

midacoRprint <- function( C,printeval,save2file,iflag,istop,f,g,x,xl,xu,
                          o,n,ni,m,me,rw,pf,maxeval,maxtime,param,P,T,key)
{
########################################################################
########################################################################
#######################################################################
#     This subroutine handles all printing commands for MIDACO.
#     This subroutine will initialize IFLAG=0 and ISTOP=0.
#     This subroutine will check the MAXEVAL and MAXTIME criteria.  
########################################################################
########################################################################
########################################################################    
 if( C == 1 )
 {
  #################################    
  bestg <<- double(length = m+2*o)
  bestx <<- double(length = n)  
  #################################
  iflag = 0
  istop = 0
  tmax <<- as.double(maxtime)  
  # tstart <<- Sys.time()  
  timearray = proc.time()
  tstart <<- timearray[3]

  Reval <<- 0
  if( param[1] <= 0.0 )
  {
    acc <<- 1.0E-3
  }else{
    acc <<- param[1]
  }
  EXTRAOFFSET <<- 5*(P+o+n+m)+100
  QQQ  <<- 102*n+(m+2*o)+516 + EXTRAOFFSET
  kx   <<- 9          
  kf   <<- 9+1+n
  kg   <<- 9+1+n
  kres <<- 9+1+n+1+ m
  if( o > 1 ){ kres <<- 9+1+n+1+ (m+2*o) }
  wx   <<- QQQ
  wf   <<- QQQ+1+n
  wg   <<- QQQ+1+n
  wres <<- QQQ+1+n+1+ m
  if( o > 1 ){ wres <<- QQQ+1+n+1+ (m+2*o) }
  if( (save2file >= 1) && (printeval >= 1) )
  {
    fn <- "MIDACO_SCREEN.TXT"
    if (file.exists(fn)) file.remove(fn)
    fn <- "MIDACO_SOLUTION.TXT"
    if (file.exists(fn)) file.remove(fn)    
  }
  bestf <<- 1.0E+99
  bestr <<- 1.0E+99
  dummy_f   <<- 1.0E+99
  dummy_vio <<- 1.0E+99   
  tic <<- 0
  pfmax <<- 1000
  if(param[10] >=  1.0){ pfmax <<-  as.integer( param[10]) }
  if(param[10] <= -1.0){ pfmax <<-  as.integer(-param[10]) }
  if(printeval >= 1)
  {
      print_head( P, o, n,ni,m,me, param, maxeval, maxtime,
                  printeval, save2file, key, 0)
      if(save2file >= 1)
      {
       print_head( P, o, n,ni,m,me, param, maxeval, maxtime,
                   printeval, save2file, key, 1)        
      }
      if(save2file >= 1)
      {
ABC <- sprintf("
MIDACO - SOLUTION
-----------------
This file saves the current best solution X found by MIDACO.
This file is updated after every PRINTEVAL function evaluation,
if X has been improved.\n \n")        
      }   
      cat(ABC, file="MIDACO_SOLUTION.TXT", append=TRUE) 
  }
    
 flush.console()
 }
########################################################################
########################################################################
########################################################################
 if( C == 2 )
 {
    # tnow = Sys.time() - tstart     
    timearray = proc.time()
    tnow = timearray[3] - tstart
   
    Reval <<- Reval + P
    if( iflag >= 10 )
    {
        warnings_and_erros( iflag, 0 )
        if( save2file >= 1 ){ warnings_and_erros( iflag, 1 ) }
        if( save2file >= 1 ){ warnings_and_erros( iflag, 2 ) }
    }
    if( printeval >= 1 )
    {
        tic <<- tic + P
        if( tic >= printeval || Reval == P || iflag >= 1)
        {
            if( Reval < 0){
                tic <<- 0
                Reval <<- -Reval-2*printeval
            }
            if( Reval > P ){ tic <<- 0 }
            if( rw[kres] == rw[wres] )
            {
                kbest = kf
                wbest = wf
            }
            else
            {
                kbest = kres
                wbest = wres
            }
            if( rw[wbest] < rw[kbest] || 
                (iflag>=1||iflag==-300||iflag==-500))
            {
                bestf = rw[wf]
                bestr = rw[wres]
                if( o <= 1 )
                {
                    for (i in 1:m){ bestg[i] = rw[wg+i] }
                }
                else
                {
                    for (i in 1:(m+2*o)){ bestg[i] = rw[wg+i] }                 
                }
                for (i in 1:n){ bestx[i] = rw[wx+i] }
            }
            else
            {
                bestf = rw[kf]
                bestr = rw[kres]
                if( o <= 1 )
                {
                    for (i in 1:m){ bestg[i] = rw[kg+i] }
                }
                else
                {
                    for (i in 1:(m+2*o)){ bestg[i] = rw[kg+i] }    
                }
                for (i in 1:n){ bestx[i] = rw[kx+i] }                
            }
            psize = as.integer(pf[1])
            if( iflag < 100 )
            {
                print_line(o,Reval,tnow,bestf,bestr,psize,0) 
                if( save2file >= 1 ){ print_line(o,Reval,tnow,bestf,bestr,psize,1) }
            }
            if( save2file >= 1 )
            {
                update = 0
                if( (bestr<dummy_vio) || (bestr==dummy_vio)&&(bestf<dummy_f-1.0e-12) )
                {
                    dummy_f <<- bestf
                    dummy_vio <<- bestr
                    update = 1
                }
                if( update == 1 )
                {
                    ABC <- sprintf("\n\n            CURRENT BEST SOLUTION") 
                    cat(ABC, file="MIDACO_SOLUTION.TXT", append=TRUE)   
                    print_solution( o,n,m,me, bestx,bestg,bestf,bestr, pf, xl,xu,acc,Reval,tnow, iflag, 2)                                  
                }

            }
            if( save2file >= 1 )
            {
                if( o > 1 ){ print_paretofront(o,m,n,pf,bestx,bestg,pfmax,iout) }
            }
        }
        if( istop == 0 )
        {
            if( tnow >= tmax ){ iflag = -999 }
            if( Reval >= maxeval-1){ if( maxeval < 999999999 ){ iflag = -999 } }
        }
        if( istop >= 1 )
        {

            bestf = f[1]
            if( o >= 2 ){ bestf = rw[wf] }

            for (i in 1:m){ bestg[i] = g[i] }
            if( o >= 2 )
            {
                for (i in 1:(o*2)){ bestg[m+i] = rw[wg+m+i] }
            }
            if( printeval > 0 )
            {
        

                print_final(iflag,tnow,tmax,Reval,maxeval,o,n,m,me,
                            x,bestg,bestf,xl,xu,rw,acc,wres,pf,param, 0 )
                if( save2file >= 1 )
                {
                    print_final(iflag,tnow,tmax,Reval,maxeval,o,n,m,me,
                                x,bestg,bestf,xl,xu,rw,acc,wres,pf,param, 1 )
                    print_final(iflag,tnow,tmax,Reval,maxeval,o,n,m,me,
                                x,bestg,bestf,xl,xu,rw,acc,wres,pf,param, 2 )                    
                }                
            }
        }
    }
  flush.console()
 } 
########################################################################
########################################################################
########################################################################
  # dummy = 0
  output <- list( iflag )
  return(output) 
}




########################################################################
########################################################################
########################################################################
print_head <- function( P, o, n,ni,m,me, par, maxeval, maxtime,
                        printeval, save2file, key, iout){

ABC <- sprintf("
 MIDACO 6.0    (www.midaco-solver.com)
 -------------------------------------\n
 LICENSE-KEY:  %60s\n
 ----------------------------------------
 | OBJECTIVES%5i | PARALLEL%10i |
 |--------------------------------------|
 | N%14i | MAXEVAL%11i |
 | NI%13i | MAXTIME%11i |     
 | M%14i | PRINTEVAL%9i |     
 | ME%13i | SAVE2FILE%9i |
 |--------------------------------------|",
 key,o,P,n,maxeval,ni,maxtime,m,printeval,me,save2file)

 if( iout == 0 ){ cat(ABC) }
 if( iout == 1 ){ cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE) }

 d=0

 if( iout == 0 )
 {
  if( par[ 1] > 0.0 | par[ 1] < 0.0 ){ ABC <- sprintf("\n | PARAM( 1) %14.7e ACCURACY    |",par[ 1]); cat(ABC); d=1 }
  if( par[ 2] > 0.0 | par[ 2] < 0.0 ){ ABC <- sprintf("\n | PARAM( 2) %14.7e RANDOM-SEED |",par[ 2]); cat(ABC); d=1 }
  if( par[ 3] > 0.0 | par[ 3] < 0.0 ){ ABC <- sprintf("\n | PARAM( 3) %14.7e FSTOP       |",par[ 3]); cat(ABC); d=1 }
  if( par[ 4] > 0.0 | par[ 4] < 0.0 ){ ABC <- sprintf("\n | PARAM( 4) %14.7e ALGOSTOP    |",par[ 4]); cat(ABC); d=1 }
  if( par[ 5] > 0.0 | par[ 5] < 0.0 ){ ABC <- sprintf("\n | PARAM( 5) %14.7e EVALSTOP    |",par[ 5]); cat(ABC); d=1 }
  if( par[ 6] > 0.0 | par[ 6] < 0.0 ){ ABC <- sprintf("\n | PARAM( 6) %14.7e FOCUS       |",par[ 6]); cat(ABC); d=1 }
  if( par[ 7] > 0.0 | par[ 7] < 0.0 ){ ABC <- sprintf("\n | PARAM( 7) %14.7e ANTS        |",par[ 7]); cat(ABC); d=1 }
  if( par[ 8] > 0.0 | par[ 8] < 0.0 ){ ABC <- sprintf("\n | PARAM( 8) %14.7e KERNEL      |",par[ 8]); cat(ABC); d=1 }
  if( par[ 9] > 0.0 | par[ 9] < 0.0 ){ ABC <- sprintf("\n | PARAM( 9) %14.7e ORACLE      |",par[ 9]); cat(ABC); d=1 }
  if( par[10] > 0.0 | par[10] < 0.0 ){ ABC <- sprintf("\n | PARAM(10) %14.7e PARETOMAX   |",par[10]); cat(ABC); d=1 }
  if( par[11] > 0.0 | par[11] < 0.0 ){ ABC <- sprintf("\n | PARAM(11) %14.7e EPSILON     |",par[11]); cat(ABC); d=1 }
  if( par[12] > 0.0 | par[12] < 0.0 ){ ABC <- sprintf("\n | PARAM(12) %14.7e BALANCE     |",par[12]); cat(ABC); d=1 }
  if( par[13] > 0.0 | par[13] < 0.0 ){ ABC <- sprintf("\n | PARAM(13) %14.7e CHARACTER   |",par[13]); cat(ABC); d=1 }  
 }
 if( iout == 1 )
 {
  if( par[ 1] > 0.0 | par[ 1] < 0.0 ){ ABC <- sprintf("\n | PARAM( 1) %14.7e ACCURACY    |",par[ 1]); cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE); d=1 }
  if( par[ 2] > 0.0 | par[ 2] < 0.0 ){ ABC <- sprintf("\n | PARAM( 2) %14.7e RANDOM-SEED |",par[ 2]); cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE); d=1 }
  if( par[ 3] > 0.0 | par[ 3] < 0.0 ){ ABC <- sprintf("\n | PARAM( 3) %14.7e FSTOP       |",par[ 3]); cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE); d=1 }
  if( par[ 4] > 0.0 | par[ 4] < 0.0 ){ ABC <- sprintf("\n | PARAM( 4) %14.7e ALGOSTOP    |",par[ 4]); cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE); d=1 }
  if( par[ 5] > 0.0 | par[ 5] < 0.0 ){ ABC <- sprintf("\n | PARAM( 5) %14.7e EVALSTOP    |",par[ 5]); cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE); d=1 }
  if( par[ 6] > 0.0 | par[ 6] < 0.0 ){ ABC <- sprintf("\n | PARAM( 6) %14.7e FOCUS       |",par[ 6]); cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE); d=1 }
  if( par[ 7] > 0.0 | par[ 7] < 0.0 ){ ABC <- sprintf("\n | PARAM( 7) %14.7e ANTS        |",par[ 7]); cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE); d=1 }
  if( par[ 8] > 0.0 | par[ 8] < 0.0 ){ ABC <- sprintf("\n | PARAM( 8) %14.7e KERNEL      |",par[ 8]); cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE); d=1 }
  if( par[ 9] > 0.0 | par[ 9] < 0.0 ){ ABC <- sprintf("\n | PARAM( 9) %14.7e ORACLE      |",par[ 9]); cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE); d=1 }
  if( par[10] > 0.0 | par[10] < 0.0 ){ ABC <- sprintf("\n | PARAM(10) %14.7e PARETOMAX   |",par[10]); cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE); d=1 }
  if( par[11] > 0.0 | par[11] < 0.0 ){ ABC <- sprintf("\n | PARAM(11) %14.7e EPSILON     |",par[11]); cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE); d=1 }
  if( par[12] > 0.0 | par[12] < 0.0 ){ ABC <- sprintf("\n | PARAM(12) %14.7e BALANCE     |",par[12]); cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE); d=1 }
  if( par[13] > 0.0 | par[13] < 0.0 ){ ABC <- sprintf("\n | PARAM(12) %14.7e CHARACTER   |",par[13]); cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE); d=1 }
 } 

 if( d == 0 )
 {
    ABC <- sprintf("\n | PARAMETER:    All by default (0)     |") 
    if( iout == 0 ){ cat(ABC) }
    if( iout == 1 ){ cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE) } 
 }

 ABC <- sprintf("\n ----------------------------------------")
 if( iout == 0 ){ cat(ABC) }
 if( iout == 1 ){ cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE) } 


if( o == 1 )
{
ABC <- sprintf(
"\n \n
 [     EVAL,    TIME]        OBJECTIVE FUNCTION VALUE         VIOLATION OF G(X)
 ------------------------------------------------------------------------------
")
if( iout == 0 ){ cat(ABC) }
if( iout == 1 ){ cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE) } 
}
else
{
ABC <- sprintf(
"\n \n
 [     EVAL,    TIME]   MULTI-OBJECTIVE PROGRESS   VIOLATION OF G(X)
 -------------------------------------------------------------------   [PARETO]
")
if( iout == 0 ){ cat(ABC) }
if( iout == 1 ){ cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE) } 
}

}
########################################################################
########################################################################
########################################################################
print_line <- function(o,Reval,tnow,f,vio,psize,iout)
{
    if( o == 1 )
    {
        if( abs(f) <= 1.0E+10 )
        {
            if( vio <= 1.0E+5 )
            {
              ABC <- sprintf(' [%9i,%8.2f]        F(X):%19.8f         VIO:%13.6f\n',Reval,tnow,f,vio)
            }
            else
            {
              ABC <- sprintf(' [%9i,%8.2f]        F(X):%19.8f         VIO:%13.6e\n',Reval,tnow,f,vio)
            }
        }
        else
        {
            if( vio <= 1.0E+5 )
            {
              ABC <- sprintf(' [%9i,%8.2f]        F(X):%19.8e         VIO:%13.6f\n',Reval,tnow,f,vio)
            }
            else
            {
              ABC <- sprintf(' [%9i,%8.2f]        F(X):%19.8e         VIO:%13.6e\n',Reval,tnow,f,vio)
            }            
        }
    }
    else
    {
        if( abs(f) <= 1.0E+10 )
        {
            if( vio <= 1.0E+5 )
            {
              ABC <- sprintf(' [%9i,%8.2f]   PRO:%20.8f   VIO:%13.6f   [%6i]\n',Reval,tnow,f,vio,psize)
            }
            else
            {
              ABC <- sprintf(' [%9i,%8.2f]   PRO:%20.8f   VIO:%13.6e   [%6i]\n',Reval,tnow,f,vio,psize)
            }
        }
        else
        {
            if( vio <= 1.0E+5 )
            {
              ABC <- sprintf(' [%9i,%8.2f]   PRO:%20.8e   VIO:%13.6f   [%6i]\n',Reval,tnow,f,vio,psize)
            }
            else
            {
              ABC <- sprintf(' [%9i,%8.2f]   PRO:%20.8e   VIO:%13.6e   [%6i]\n',Reval,tnow,f,vio,psize)
            }            
        }
    }
    ####################################################################
    ####################################################################
    ####################################################################
    if( iout == 0 ){ cat(ABC) }
    if( iout == 1 ){ cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE) }

}
########################################################################
########################################################################
########################################################################
print_solution <- function(o,n,m,me, x,g,f,vio, pf, xl,xu,acc,Reval,tnow, iflag, iout)
{
ABC <- sprintf("
 --------------------------------------------
 EVAL:%9i,  TIME:%9.2f,  IFLAG:%4i
 --------------------------------------------",Reval,tnow,iflag)

 if( iout == 0 ){ cat(ABC) }
 if( iout == 1 ){ cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE) }
 if( iout == 2 ){ cat(ABC, file="MIDACO_SOLUTION.TXT", append=TRUE) }

if( o == 1 )
{
ABC <- sprintf(" 
 F(X) =%38.15g
",f)
if( iout == 0 ){ cat(ABC) }
if( iout == 1 ){ cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE) }
if( iout == 2 ){ cat(ABC, file="MIDACO_SOLUTION.TXT", append=TRUE) }
}
else
{

psize = as.integer(pf[1])

ABC <- sprintf("
 PROGRESS%36.15g
 --------------------------------------------
 NUMBER OF PARETO POINTS%21i
 --------------------------------------------
",f,psize) 

  if( iout == 0 ){ cat(ABC) }
  if( iout == 1 ){ cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE) }
  if( iout == 2 ){ cat(ABC, file="MIDACO_SOLUTION.TXT", append=TRUE) }
 for (i in 1:o)
 {
    if( g[m+o+i] <= 0.0 ){ objvalue =  g[m+i] }
    if( g[m+o+i]  > 0.0 ){ objvalue = -g[m+i] }
    ABC <- sprintf(" F(%3i) = %35.15g\n",i,objvalue)
    if( iout == 0 ){ cat(ABC) }
    if( iout == 1 ){ cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE) }
    if( iout == 2 ){ cat(ABC, file="MIDACO_SOLUTION.TXT", append=TRUE) }
 }  

} # end of else clause

if( m > 0 ){

ABC <- sprintf(" --------------------------------------------
 VIOLATION OF G(X) %26.12g
 -------------------------------------------- 
",vio)
if( iout == 0 ){ cat(ABC) }
if( iout == 1 ){ cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE) }
if( iout == 2 ){ cat(ABC, file="MIDACO_SOLUTION.TXT", append=TRUE) }  
for (i in 1:m)
{
if( m <= me )
{ 
    if( abs(g[i]) <= acc){ ABC <- sprintf(" G(%4i) =%16.8g  (EQUALITY CONSTR)\n",i,g[i]) }
    if( abs(g[i])  > acc){ ABC <- sprintf(" G(%4i) =%16.8g  (EQUALITY CONSTR)  <---  INFEASIBLE  ( G NOT = 0 )\n",i,g[i]) }    
}
if( m  > me )
{ 
    if( g[i] >= -acc){ ABC <- sprintf(" G(%4i) =%16.8g  (IN-EQUAL CONSTR)\n",i,g[i]) }
    if( g[i]  < -acc){ ABC <- sprintf(" G(%4i) =%16.8g  (IN-EQUAL CONSTR)  <---  INFEASIBLE  ( G < 0 )\n",i,g[i]) }    
}
 if( iout == 0 ){ cat(ABC) }
 if( iout == 1 ){ cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE) }
 if( iout == 2 ){ cat(ABC, file="MIDACO_SOLUTION.TXT", append=TRUE) }
}

}




ABC <- sprintf(" --------------------------------------------            BOUNDS-PROFIL\n") 
 if( iout == 0 ){ cat(ABC) }
 if( iout == 1 ){ cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE) }
 if( iout == 2 ){ cat(ABC, file="MIDACO_SOLUTION.TXT", append=TRUE) }


for (i in 1:n)
{

profil = -1

if( x[i] > xu[i] + 1.0e+6 ){ profil = 91 }
if( x[i] < xl[i] - 1.0e+6 ){ profil = 92 }
if( xl[i] > xu[i] ){ profil = 93 }
if( abs(xl[i]-xu[i]) < 1.0e+8 ){ profil = 90 }

breakloop = 0
for (k in 1:21)
{
    if( breakloop == 0 ){
    if( x[i] <= xl[i] + as.double(k) * (xu[i]-xl[i])/21.0 )
    {
        profil = k
        breakloop = 1
    }}
}

if( abs(x[i]-xl[i]) <= (xu[i]-xl[i])/1000.0 ){ profil = 0}
if( abs(x[i]-xu[i]) <= (xu[i]-xl[i])/1000.0 ){ profil = 22}

if( profil ==  0 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   XL___________________\n",i,x[i]) }
if( profil ==  1 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   x____________________\n",i,x[i]) }
if( profil ==  2 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   _x___________________\n",i,x[i]) }
if( profil ==  3 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   __x__________________\n",i,x[i]) }
if( profil ==  4 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   ___x_________________\n",i,x[i]) }
if( profil ==  5 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   ____x________________\n",i,x[i]) }
if( profil ==  6 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   _____x_______________\n",i,x[i]) }
if( profil ==  7 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   ______x______________\n",i,x[i]) }
if( profil ==  8 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   _______x_____________\n",i,x[i]) }
if( profil ==  9 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   ________x____________\n",i,x[i]) }
if( profil == 10 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   _________x___________\n",i,x[i]) }
if( profil == 11 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   __________x__________\n",i,x[i]) }
if( profil == 12 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   ___________x_________\n",i,x[i]) }
if( profil == 13 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   ____________x________\n",i,x[i]) }
if( profil == 14 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   _____________x_______\n",i,x[i]) }
if( profil == 15 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   ______________x______\n",i,x[i]) }
if( profil == 16 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   _______________x_____\n",i,x[i]) }
if( profil == 17 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   ________________x____\n",i,x[i]) }
if( profil == 18 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   _________________x___\n",i,x[i]) }
if( profil == 19 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   __________________x__\n",i,x[i]) }
if( profil == 20 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   ___________________x_\n",i,x[i]) }
if( profil == 21 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   ____________________x\n",i,x[i]) }
if( profil == 22 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   ___________________XU\n",i,x[i]) }
if( profil == 90 ){ ABC <- sprintf(" X(%5i) = %33.16g    #   WARNING: XL = XU     \n",i,x[i]) }
if( profil == 91 ){ ABC <- sprintf(" X(%5i) = %33.16g  <---  *** ERROR *** (X > XU) \n",i,x[i]) }
if( profil == 92 ){ ABC <- sprintf(" X(%5i) = %33.16g  <---  *** ERROR *** (X < XL) \n",i,x[i]) }
if( profil == 93 ){ ABC <- sprintf(" X(%5i) = %33.16g  <---  *** ERROR *** (XL > XU)\n",i,x[i]) }

 if( iout == 0 ){ cat(ABC) }
 if( iout == 1 ){ cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE) }
 if( iout == 2 ){ cat(ABC, file="MIDACO_SOLUTION.TXT", append=TRUE) }

} # end of forloop

ABC <- sprintf("\n")

 if( iout == 0 ){ cat(ABC) }
 if( iout == 1 ){ cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE) }
 if( iout == 2 ){ cat(ABC, file="MIDACO_SOLUTION.TXT", append=TRUE) } 

}    
########################################################################
########################################################################
########################################################################
print_final <- function(iflag,tnow,tmax,eval,maxeval,o,n,m,me,
                        x,g,f,xl,xu,rw,acc,wres,pf,param, iout )
{

    if( iflag == 1 || iflag == 2 )
    {
        if( tnow >= tmax ){     ABC <- sprintf("\n OPTIMIZATION FINISHED  --->  MAXTIME REACHED") }
        if( Reval >= maxeval ){ ABC <- sprintf("\n OPTIMIZATION FINISHED  --->  MAXEVAL REACHED") }            
    }
    if( iflag == 3 || iflag == 4 )
    {
        ABC <- sprintf("\n OPTIMIZATION FINISHED  --->  ALGOSTOP (=%3i)",param[4])           
    }  
    if( iflag == 5 || iflag == 6 )
    {
        ABC <- sprintf("\n OPTIMIZATION FINISHED  --->  EVALSTOP (=%9i)",param[5])           
    }  
    if( iflag == 7 )
    {
        ABC <- sprintf("\n OPTIMIZATION FINISHED  --->  FSTOP REACHED")           
    }            


    if( iflag <= 10 )
    {
      if( iout == 0 ){ cat(ABC) }
      if( iout == 1 ){ cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE) }
      if( iout == 2 ){ cat(ABC, file="MIDACO_SOLUTION.TXT", append=TRUE) }
    }    

    ABC <- sprintf("\n\n\n         BEST SOLUTION FOUND BY MIDACO       ")

    if( iout == 0 ){ cat(ABC) }
    if( iout == 1 ){ cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE) }
    if( iout == 2 ){ cat(ABC, file="MIDACO_SOLUTION.TXT", append=TRUE) }

    print_solution( o,n,m,me, x,g,f,rw[wres], pf, xl,xu,acc,Reval,tnow, iflag, iout)
}




warnings_and_erros <- function( iflag, iout)
{
    if( iflag < 100 )
    {
        ABC <- sprintf("\n *** WARNING ***   ( IFLAG = %i ) \n",iflag) 
    }
    else
    {
        ABC <- sprintf("\n\n\n *** *** MIDACO INPUT ERROR *** ***   ( IFLAG = %i ) \n\n\n",iflag)
    }
    if( iout == 0 ){ cat(ABC) }
    if( iout == 1 ){ cat(ABC, file="MIDACO_SCREEN.TXT", append=TRUE) }
    if( iout == 2 ){ cat(ABC, file="MIDACO_SOLUTION.TXT", append=TRUE) }        
}




print_paretofront <- function(o,m,n,pf,bestx,bestg,pfmax,iout)
{
ABC <- sprintf("#########################################################
### This file contains the pareto front approximation ###
#########################################################
### Solution format:     F(1:O)    G(1:M)    X(1:N)   ###
#########################################################\n")

cat(ABC, file="MIDACO_PARETOFRONT.tmp", append=FALSE)

psize = as.integer(pf[1])

ABC <- sprintf("#
#        O         M         N     PSIZE
#
 %9i %9i %9i %9i\n",o,m,n,psize)    

cat(ABC, file="MIDACO_PARETOFRONT.tmp", append=TRUE)

ABC <- sprintf("#
#        MIDACO solution
#\n")    

cat(ABC, file="MIDACO_PARETOFRONT.tmp", append=TRUE)

    for (i in 1:o)
    {
        dummy = bestg[m+i]
        if( bestg[m+o+i] > 0.0 ){ dummy = -dummy }
        ABC <- sprintf(" %18.7g",dummy)
        cat(ABC, file="MIDACO_PARETOFRONT.tmp", append=TRUE)
    }
    if( m > 0 ){
    for (i in 1:m)
    {
        dummy = bestg[i]
        ABC <- sprintf(" %18.7g",dummy)
        cat(ABC, file="MIDACO_PARETOFRONT.tmp", append=TRUE)
    }
    }   
    for (i in 1:n)
    {
        dummy = bestx[i]
        ABC <- sprintf(" %18.7g",dummy)
        cat(ABC, file="MIDACO_PARETOFRONT.tmp", append=TRUE)
    }        

    ABC <- sprintf("\n")
    cat(ABC, file="MIDACO_PARETOFRONT.tmp", append=TRUE)

ABC <- sprintf("#
#        All non-dominated solutions found by MIDACO
#\n")    

cat(ABC, file="MIDACO_PARETOFRONT.tmp", append=TRUE)

if( psize > 0 ){
for (k in 1:psize)
{
    for (i in 1:o)
    {
        dummy = pf[ 2 + o*(k-1)+i-1 ]
        ABC <- sprintf(" %18.7g",dummy)
        cat(ABC, file="MIDACO_PARETOFRONT.tmp", append=TRUE)
    }
    if( m > 0 ){
    for (i in 1:m)
    {
        dummy = pf[ 2 + o*pfmax + m*(k-1)+i-1 ]
        ABC <- sprintf(" %18.7g",dummy)
        cat(ABC, file="MIDACO_PARETOFRONT.tmp", append=TRUE)
    }   
    }
    for (i in 1:n)
    {
        dummy = pf[ 2 + o*pfmax + m*pfmax + n*(k-1)+i-1 ]
        ABC <- sprintf(" %18.7g",dummy)
        cat(ABC, file="MIDACO_PARETOFRONT.tmp", append=TRUE)
    }        

    ABC <- sprintf("\n")
    cat(ABC, file="MIDACO_PARETOFRONT.tmp", append=TRUE)
}
}

file.rename( "MIDACO_PARETOFRONT.tmp", "MIDACO_PARETOFRONT.TXT" )

}

########################################################################
############################ END OF FILE ###############################
########################################################################