adjust.y = function(df,yform,comps,type="FIXED"){
  
  #' step1: create outputs
  ngid=nlevels(df$gid)
  yform=as.formula(yform)
  
  #' step2: provide computations
  if(type=="FIXED"){m = try(lm(yform,df))}
  if(type=="MIXED"){m = try(lmer(yform,df))}
  
  #' step3: returning output
  y = data.frame(try(emmeans(m,specs = comps)))
  return(y)
}

## modelo
#form = "yield~gid*Env+(1|Stand:gid)+(1|MACRO)"
#form = "yield~(1|gid)+Env+(1|gid:Env)+stand"
#ncomps=c("gid","Env")
#y=adjust.y(df = data,yform = form,comps = ncomps,type = "MIXED")


pheno.rep = function(Gid,Env,value,df){
  pheno<-c()
  
  data = df[,c(Gid,Env,value)]
  names(data)= c("Gid","Env","value")
  for(i in 1:nlevels(data$Gid)){
    
    data2<-droplevels(data[-which(data$Gid == levels(data$Gid)[i]),])
    
    for(j in 1:nlevels(data2$Env)){
      data3<-droplevels(data2[-which(data2$Env == levels(data2$Env)[j]),])
      
      m.prod<-lmer(value~(1|Gid)+(1|Env:Gid),data3)
      
      ge<-data.frame(ranef(m.prod)[1])                # interacao genotipo x local
      
      g<-data.frame(ranef(m.prod)[2])
      g$g<-row.names(g)
      names(g)[1]<-"g.hat"
      
      GE<-data.frame(colsplit(row.names(ge), ":",c("e","g")),ge.hat=ge$X.Intercept.)
      names(GE)
      
      Y<-merge(g,GE,by="g")
      head(Y)
      Y$value.trait <- Y$g.hat+Y$ge.hat+fixef(m.prod)[1]
      pheno<-rbind(pheno,Y)
    }
  } 
  return(pheno)
}


FW.out = function(y, gid, env, method = "OLS",scale=TRUE, to = c(0.5,2)){
  fit = FW(y=y,VAR = gid,ENV = env,method = method)
  out = data.frame(g = fit$mu+fit$g, b = fit$b+1)
  if(scale == TRUE){out$b = rescale(out$b,to=to)}
  names(out) = c("g","b")
  out$gid = rownames(out)
  return(list(acc=cor(fit$y,fit$yhat),adp=out))
}


ge.means = function(k){
  ge = as.data.frame(ranef(k)[1])
  ge = data.frame(ge,colsplit(rownames(ge),pattern = ":",c("gid","Env")))
  names(ge)[1] = "ge"
  
  e = as.data.frame(ranef(k)[2])
  e$Env = rownames(e)
  names(e)[1] = "e"
  
  g = data.frame(emmeans(k,specs = c("gid")))
  
  y = merge(g,merge(ge,e,by="Env"),by="gid")
  y$yield = y$emmean + y$e + y$ge
  return(y)
}

