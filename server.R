require(shiny)
require(fishmethods)
require(ggplot2)
require(truncnorm)
require(data.table)
require(RColorBrewer) 
require(viridis)

#shinyServer(
  function(input, output) 
  {    
    
    Then_M<-function(Amax)
    {
      M_vals<-rep(NA,4)
      M_vals[1]<-4.889*Amax^-0.916
      M_vals[2]<-5.109/Amax
      M_vals[3]<-exp(1.717-1.01*log(Amax))
      M_vals[4]<-5.4/Amax
      return(M_vals)
    }
    
    Then_VBGF<-function(Linf,k)
    {
      M_val_vbgf<-4.11*k^0.73*Linf^-0.33
      return(M_val_vbgf)
    }
    
    Jensen_M_amat<-function(Amat)
    {
      M_val_Jensen<-1.65/Amat
      return(M_val_Jensen)
    }
    
    Jensen_M_k<-function(k)
    {
      M_val_Jensen_k<-k*c(1.5,1.6)
      return(M_val_Jensen_k)
    }
    
    #Rikhter & Efanov
    Rikhter_Efanov_Amat_M<-function(Amat)
      {
        M_val_RiEf<-(1.52/Amat^0.72)-0.16
        return(M_val_RiEf)
      }
    
        
    Chen_N_Wat_M<-function(Amax,k,t0,out.type=1)
    {
      if(anyNA(c(Amax,k,t0))){M.out<-NA}
      else
      {
        M_ages<-rep(NA,length(c(1:Amax)))
        tM<--1/k*(log(abs(1-exp(k*t0)))+t0)
        a0<-1-exp(-k*(tM-t0))
        a1<-k*exp(-k*(tM-t0))
        a2<--0.5*k^2*exp(-k*(tM-t0))
        for(a in 1:Amax)
        {
          if(a<=tM){M_ages[a]<-k/(1-exp(-k*(a-t0)))}
          if(a>tM){M_ages[a]<-k/(a0+a1*(a-tM)+a2*(a-tM)^2)}
        }
        if(out.type==1){M.out<-mean(M_ages)}
        else M.out<-M_ages
      }
      return(M.out)
    }

    Gislason_M_a<-function(Amax,Linf,k,t0)
    {
		Lts<-Linf*(1-exp(-k*(c(1:Amax)-t0)))
		Gis_Ms_a<-mapply(function(x) M.empirical(Linf=Linf,Kl=k,Bl=Lts[x],method=9)[1],x=1:length(Lts),SIMPLIFY=TRUE)
		return(Gis_Ms_a)
    }

	fishlife.M <- function(species){
 	 # Setup container
 	 if(species!="")
 	 {
    	 spp <- sort(unique(species))
     	fl <- data.frame(species=spp, linf_cm=NA, k=NA, winf_g=NA, tmax_yr=NA, tmat_yr=NA,
                      m=NA, lmat_cm=NA, temp_c=NA, stringsAsFactors=F)
     
   	 # Loop through species
     	for(i in 1:nrow(fl)){
        
        # Get spp info
        sciname <- fl$species[i]
        genus <- stringr::word(sciname, 1)
        nwords_in_spp <- length(strsplit(sciname, " ")[[1]])
        species <- stringr::word(sciname, start=2, end=nwords_in_spp)
        species <- ifelse(species=="spp", "predictive", species)
        
        # Try looking up in FishLife
        spp_info <- try(FishLife::Plot_taxa(FishLife::Search_species(Genus=genus, Species=species)$match_taxonomy,mfrow=c(2,2)))
        if(inherits(spp_info, "try-error")){
           # Record blanks
        }else{
           # Values are in log-scale except temperature
           spp_lh_vals_log <- spp_info[[1]]$Mean_pred
           spp_lh_vals <- c(exp(spp_lh_vals_log[1:7]), spp_lh_vals_log[8],spp_lh_vals_log[9:20])
       	 }
     	}
     
     # Return
     	return(as.numeric(spp_lh_vals)[6])
  	}
  	 else{return(NA)}     
	}
    
####### END FUNCTIONS ########
    
    
 M_vals_all<- reactive({
   fishlife.M.out<-Pauly80lt_M<-Pauly80wt_M<-AnC75_M<-Roff_M<-GnD_GSI_M<-PnW_M<-Lorenzen96_M<-Gislason_M<-Gislason_M_ages<-NA
   
   if(input$Genspp!="Type Genus and species here"){fishlife.M.out<-fishlife.M(input$Genspp)}
   Then_M_Amax<-Then_M(input$Amax)
   if(!(anyNA(c(input$k_vbgf,input$Amax)))){AnC75_M<-M.empirical(Kl=input$k_vbgf,tmax=input$Amax,method=4)[1]}
   Then_M_VBGF<-Then_VBGF(input$Linf*10,input$k_vbgf)
   Jensen_M_VBGF<-Jensen_M_k(input$k_vbgf) 
   if(!(anyNA(c(input$Linf,input$k_vbgf,input$Lt_in))))
   	{
   		Gislason_M<-M.empirical(Linf=input$Linf,Kl=input$k_vbgf,Bl=input$Lt_in,method=9)[1]
   	}
   if(!(anyNA(c(input$Amax,input$Linf,input$k_vbgf,input$Lt_in,input$t0))))
   	{
   		Gislason_M_ages<-Gislason_M_a(input$Amax,input$Linf,input$k_vbgf,input$t0)
   	}
   CnW_M_VBGF<-Chen_N_Wat_M(input$Age_in,input$k_vbgf,input$t0)
   CnW_M_a_VBGF<-Chen_N_Wat_M(input$Amax,input$k_vbgf,input$t0,out.type = 0)
   maxage<-input$Amax
   if(!is.na(maxage)){CnW_M_a_VBGF_table<-cbind(c(1:maxage),CnW_M_a_VBGF,Gislason_M_ages)
   colnames(CnW_M_a_VBGF_table)<-c("Age","CnW_M","Gislason_M")}
   if(!(anyNA(c(input$k_vbgf,input$Amat)))){Roff_M<-M.empirical(Kl=input$k_vbgf,tm=input$Amat,method=5)[1]}
   Jensen_M_Amat<-Jensen_M_amat(input$Amat)
   Rikhter_Efanov_Amat<-Rikhter_Efanov_Amat_M(input$Amat)
   if(!(anyNA(c(input$Wdry)))){PnW_M<-M.empirical(Wdry=input$Wdry,method=7)[1]}
   if(!(anyNA(c(input$Wwet)))){Lorenzen96_M<-M.empirical(Wwet=input$Wwet,method=8)[1]}
   if(!(anyNA(c(input$Linf,input$k_vbgf,input$Temp)))){Pauly80lt_M<-M.empirical(Linf=input$Linf,Kl=input$k_vbgf,TC=input$Temp,method=1)[1]}
   if(!(anyNA(c(input$Winf,input$kw,input$Temp)))){Pauly80wt_M<-M.empirical(Winf=input$Winf,Kw=input$kw,TC=input$Temp,method=2)[1]}
   if(!(anyNA(c(input$GSI)))){GnD_GSI_M<-M.empirical(GSI=input$GSI,method=6)[1]}
   User_M<-as.numeric(trimws(unlist(strsplit(input$User_M,","))))
   if(length(User_M)==0)User_M<-NA
   M_vals_all<-c(fishlife.M.out,Then_M_Amax,CnW_M_VBGF,AnC75_M,Then_M_VBGF,Jensen_M_VBGF,Gislason_M,Pauly80lt_M,Roff_M,Jensen_M_Amat,Rikhter_Efanov_Amat,Pauly80wt_M,PnW_M,Lorenzen96_M,GnD_GSI_M,User_M)
   output$downloadCW_M_a <- downloadHandler(
     filename = function() {paste0("Age_specific_M_values", '.csv') },
     content = function(file) {write.csv(CnW_M_a_VBGF_table, file=file)}
   )  
   M_vals_all
   })
     
        
   output$Mplot <- renderPlot({
   M_vals_all<-M_vals_all()
   User_M<-as.numeric(trimws(unlist(strsplit(input$User_M,","))))
   M_users<-"User input"
   if(length(User_M)>1){M_users<-paste0("User input_",c(1:length(User_M)))}
   M_methods<-c("FishLife","Then_Amax 1","Then_Amax 2","Then_Amax 3","Hamel_Amax","Chen-Wat","AnC","Then_VBGF","Jensen_VBGF 1","Jensen_VBGF 2","Gislason","Pauly_lt","Roff","Jensen_Amat","Ri_Ef_Amat","Pauly_wt","PnW","Lorenzen","GSI",M_users)
   M_types<-c("Meta-analysis",rep("Amax",4),rep("Amax:VBGF",2),rep("VBGF",4),rep("VBGF:Temp",1),"VBGF:Amat",rep("Amat",2),rep("Weight",3),rep("GSI",1),rep("User input",length(M_users)))
   M_vals_gg<-as.data.frame(cbind(M_vals_all,M_methods,M_types))
   colnames(M_vals_gg)<-c("M","Method","Input")
   M_vals_gg$Method<-factor(M_vals_gg$Method,levels=unique(M_vals_gg$Method))
   M_vals_gg$Input<-factor(M_vals_gg$Input,levels=unique(M_vals_gg$Input))     
   # plot M
   if(all(is.na(M_vals_all))){ymax<-0.5}
   if(!(all(is.na(M_vals_all)))){ymax<-ceiling((max(M_vals_all,na.rm=TRUE)*1.1*10))/10}
   
   #ggplot(M_vals_gg,aes(Method,as.numeric(Mvals),color=Mtype))+geom_point(size=2)+ylab("M")+xlab("Method")+theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))
	  if(input$M_CV==0)
	  {
	  print(ggplot(M_vals_gg,aes(Method,as.numeric(as.character(M)),color=Input))+
           geom_point(size=4)+ylab("M")+xlab("Method")+
	  	   scale_y_continuous(limits = c(0, NA))+
	  	   theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5)))
	  }
	 
      
	  if(input$M_CV>0 & input$M_CV_type=="lognormal")
	  {
      print(ggplot(M_vals_gg,aes(Method,as.numeric(as.character(M)),color=Input))+
           geom_point(size=4)+ylab("M")+xlab("Method")+
           scale_y_continuous(limits = c(0, NA))+
           theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+
   		   geom_pointrange(aes(ymin=qlnorm(0.025,log(as.numeric(as.character(M))),input$M_CV),ymax=qlnorm(0.975,log(as.numeric(as.character(M))),input$M_CV))))
      	   
      }

	if(input$M_CV>0 & input$M_CV_type=="normal")
	  {
      print(ggplot(M_vals_gg,aes(Method,as.numeric(as.character(M)),color=Input))+
           geom_point(size=4)+ylab("M")+xlab("Method")+
           theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+
   		   geom_pointrange(aes(ymin=qtruncnorm(0.025,a=0,mean=as.numeric(as.character(M)),sd=as.numeric(as.character(M))*input$M_CV),ymax=qtruncnorm(0.975,a=0,mean=as.numeric(as.character(M)),sd=as.numeric(as.character(M))*input$M_CV))))
      }
   #par(mar=c(8,4,2,6),xpd =TRUE)
   #plot(M_vals_all, col = "black",bg=c("blue","blue","blue","blue","green","green","green","green","yellow","yellow","orange","red","red","red","black","black","black","purple","brown"),xlab=" ",ylab="Natural mortality",ylim=c(0,ymax),pch=22,cex=1.5,axes=F)
   #box()
   #axis(1,at=1:length(M_vals_all),labels=M_methods,las=3)
   #axis(2)
   #legend(x="topright",legend=c("Amax","VBGF","VBGF:Temp","VBGF;Amat","Amat","Weight","GSI","User input"),pch=22,col="black",pt.bg=c("blue","green","yellow","orange","red","black","purple","brown"),bty="n",horiz=FALSE,cex=1,inset=c(-0.125,0))
   M_table<-data.frame(cbind(M_methods,M_vals_all))
   colnames(M_table)<-c("Method","M")
  # if(all(is.na(M_vals()))){return(NULL)}
   output$downloadMs <- downloadHandler(
    filename = function() {paste0("M_values", '.csv') },
    content = function(file) {write.csv(M_table, file=file)}
  )
 })

# Show the first "n" observations
 output$Mtable <- renderTable({
   fishlife.M.out<-Pauly80lt_M<-Pauly80wt_M<-AnC75_M<-Roff_M<-GnD_GSI_M<-PnW_M<-Lorenzen96_M<-Gislason_M<-NA
   if(input$Genspp!="Type Genus and species here"){fishlife.M.out<-fishlife.M(input$Genspp)}
   Then_M_Amax<-Then_M(input$Amax)
   if(!(anyNA(c(input$k_vbgf,input$Amax)))){AnC75_M<-M.empirical(Kl=input$k_vbgf,tmax=input$Amax,method=4)[1]}
   Then_M_VBGF<-Then_VBGF(input$Linf*10,input$k_vbgf)
   Jensen_M_VBGF<-Jensen_M_k(input$k_vbgf) 
   if(!(anyNA(c(input$Linf,input$k_vbgf,input$Lt_in)))){Gislason_M<-M.empirical(Linf=input$Linf,Kl=input$k_vbgf,Bl=input$Lt_in,method=9)[1]}
   CnW_M_VBGF<-Chen_N_Wat_M(input$Amax,input$k_vbgf,input$t0)
   if(!(anyNA(c(input$k_vbgf,input$Amat)))){Roff_M<-M.empirical(Kl=input$k_vbgf,tm=input$Amat,method=5)[1]}
   Jensen_M_Amat<-Jensen_M_amat(input$Amat)
   Rikhter_Efanov_Amat<-Rikhter_Efanov_Amat_M(input$Amat)
   if(!(anyNA(c(input$Wdry)))){PnW_M<-M.empirical(Wdry=input$Wdry,method=7)[1]}
   if(!(anyNA(c(input$Wwet)))){Lorenzen96_M<-M.empirical(Wwet=input$Wwet,method=8)[1]}
   if(!(anyNA(c(input$Linf,input$k_vbgf,input$Temp)))){Pauly80lt_M<-M.empirical(Linf=input$Linf,Kl=input$k_vbgf,TC=input$Temp,method=1)[1]}
   if(!(anyNA(c(input$Winf,input$kw,input$Temp)))){Pauly80wt_M<-M.empirical(Winf=input$Winf,Kw=input$kw,TC=input$Temp,method=2)[1]}
   if(!(anyNA(c(input$GSI)))){GnD_GSI_M<-M.empirical(GSI=input$GSI,method=6)[1]}
  
   M_vals_all<-c(fishlife.M.out,Then_M_Amax,CnW_M_VBGF,AnC75_M,Then_M_VBGF,Jensen_M_VBGF,Gislason_M)
   M_methods<-c("FishLife","Then_Amax 1","Then_Amax 2","Then_Amax 3","Hamel_Amax","Chen-Wat","AnC","Then_VBGF","Jensen_VBGF 1","Jensen_VBGF 2","Gislason")
   M_table<-data.frame(cbind(M_methods,signif(M_vals_all,3)))
   colnames(M_table)<-c("Method","M")
   #rownames(M_table)<-M_methods
   M_table
  })
# Show the first "n" observations
 output$Mtable2 <- renderTable({
   fishlife.M.out<-Pauly80lt_M<-Pauly80wt_M<-AnC75_M<-Roff_M<-GnD_GSI_M<-PnW_M<-Lorenzen96_M<-Gislason_M<-NA
   if(input$Genspp!="Type Genus and species here"){fishlife.M.out<-fishlife.M(input$Genspp)}
   Then_M_Amax<-Then_M(input$Amax)
   if(!(anyNA(c(input$k_vbgf,input$Amax)))){AnC75_M<-M.empirical(Kl=input$k_vbgf,tmax=input$Amax,method=4)[1]}
   Then_M_VBGF<-Then_VBGF(input$Linf*10,input$k_vbgf)
   Jensen_M_VBGF<-Jensen_M_k(input$k_vbgf) 
   if(!(anyNA(c(input$Linf,input$k_vbgf,input$Lt_in)))){Gislason_M<-M.empirical(Linf=input$Linf,Kl=input$k_vbgf,Bl=input$Lt_in,method=9)[1]}
   CnW_M_VBGF<-Chen_N_Wat_M(input$Amax,input$k_vbgf,input$t0)
   if(!(anyNA(c(input$k_vbgf,input$Amat)))){Roff_M<-M.empirical(Kl=input$k_vbgf,tm=input$Amat,method=5)[1]}
   Jensen_M_Amat<-Jensen_M_amat(input$Amat)
   Rikhter_Efanov_Amat<-Rikhter_Efanov_Amat_M(input$Amat)
   if(!(anyNA(c(input$Wdry)))){PnW_M<-M.empirical(Wdry=input$Wdry,method=7)[1]}
   if(!(anyNA(c(input$Wwet)))){Lorenzen96_M<-M.empirical(Wwet=input$Wwet,method=8)[1]}
   if(!(anyNA(c(input$Linf,input$k_vbgf,input$Temp)))){Pauly80lt_M<-M.empirical(Linf=input$Linf,Kl=input$k_vbgf,TC=input$Temp,method=1)[1]}
   if(!(anyNA(c(input$Winf,input$kw,input$Temp)))){Pauly80wt_M<-M.empirical(Winf=input$Winf,Kw=input$kw,TC=input$Temp,method=2)[1]}
   if(!(anyNA(c(input$GSI)))){GnD_GSI_M<-M.empirical(GSI=input$GSI,method=6)[1]}
   #User_M<-input$User_M
   
   M_vals_all<-c(Pauly80lt_M,Roff_M,Jensen_M_Amat,Rikhter_Efanov_Amat,Pauly80wt_M,PnW_M,Lorenzen96_M,GnD_GSI_M)
   M_methods<-c("Pauly_lt","Roff","Jensen_Amat","Ri_Ef_Amat","Pauly_wt","PnW","Lorenzen","GSI")
   M_table<-data.frame(M_vals_all)
   #rownames(M_table)<-M_methods
   #colnames(M_table)<-"M"
   M_table<-data.frame(cbind(M_methods,signif(M_vals_all,3)))
   colnames(M_table)<-c("Method","M")
   M_table
 })

 output$MtableUser <- renderTable({
   User_M<-as.numeric(trimws(unlist(strsplit(input$User_M,","))))
   if(length(User_M)==0)User_M<-NA
   print(User_M)
   M_methods<-paste0("User input_",c(1:length(User_M)))
   M_table<-data.frame(User_M)
   M_table_User<-data.frame(cbind(M_methods,signif(User_M,3)))
   colnames(M_table_User)<-c("Method","M")
   M_table_User
 })



 M.CV.method<- reactive({
  if(all(is.na(M_vals_all()))){return(NULL)}
   else{
   User_M<-as.numeric(trimws(unlist(strsplit(input$User_M,","))))
   if(length(User_M)==0)User_M<-NA
   M_users<-"User input"
   if(length(User_M)>1){M_users<-paste0("User input_",c(1:length(User_M)))}
   M.wts<-c(input$FishLife,
   			input$Then_Amax_1,
   			input$Then_Amax_2,
   			input$Then_Amax_3,
   			input$Hamel_Amax,
   			input$Chen_Wat,
   			input$AnC,
   			input$Then_VBGF,
   			input$Jensen_VBGF_1,
   			input$Jensen_VBGF_2,
   			input$Gislason,
   			input$Pauly_lt,
   			input$Roff,
   			input$Jensen_Amat,
   			input$Ri_Ef_Amat,
   			input$Pauly_wt,
   			input$PnW,
   			input$Lorenzen,
   			input$Gonosoma,
   			rep(input$UserM_wt,length(User_M)))
   names(M.wts)<-c("FishLife",
   					"Then_Amax 1",
   					"Then_Amax 2",
   					"Then_Amax 3",
   					"Hamel_Amax",
   					"Chen-Wat",
   					"AnC",
   					"Then_VBGF",
   					"Jensen_VBGF 1",
   					"Jensen_VBGF 2",
   					"Gislason",
   					"Pauly_lt",
   					"Roff",
   					"Jensen_Amat",
   					"Ri_Ef_Amat",
   					"Pauly_wt",
   					"PnW",
   					"Lorenzen",
   					"GSI",
   					M_users)
   #remove NAs
   if(any(is.na(M_vals_all()))){
     NA.ind<-attributes(na.omit(M_vals_all()))$na.action
     M.sub<-M_vals_all()[-NA.ind]
     M.wts.sub<-M.wts[-NA.ind]
   }
   else{
     M.sub<-M_vals_all()
     M.wts.sub<-M.wts
   }
   #remove 0 weight
   M.sub.n0<-M.sub[M.wts.sub>0]
   M.wts.sub.n0<-M.wts.sub[M.wts.sub>0]
   M.wts.sub.stand<-M.wts.sub.n0/sum(M.wts.sub.n0)
   samp.num<-1000000
   samps<-samp.num*M.wts.sub.stand
	
	M.CV.method <- data.table(meanval = M.sub.n0,
                       sdval = input$M_CV,
                       rep = samps,
                       wts=M.wts.sub.stand,
                       method= names(M.wts.sub.stand))
	M.CV.method
	}
	})

 	M.dists<- reactive({
		if(input$M_CV_type=="lognormal")
		{
			M.dists <- rbindlist(lapply(1:dim(M.CV.method())[1], 
                        function(x) data.table(rowval = M.CV.method()$method[x], 
                                               dist = rlnorm(M.CV.method()[x, rep],
                                               log(M.CV.method()[x, meanval]), 
                                               M.CV.method()[x, sdval]))))
		}
		if(input$M_CV_type=="normal")
		{
			M.dists <- rbindlist(lapply(1:dim(M.CV.method())[1], 
                        function(x) data.table(rowval = M.CV.method()$method[x], 
                                               dist = rtruncnorm(M.CV.method()[x, rep],
                                               a=0,
                                               b=Inf, 
                                               M.CV.method()[x, meanval], 
                                               M.CV.method()[x, meanval*sdval]))))
		}
	M.dists
 	})

#Plot Individual distributions
 output$Mdistplots<- renderPlot({ 

   if(input$M_CV==0)
   {
	 dist.dat<-M.CV.method()
	 dist.dat$method<-factor(dist.dat$method,levels=dist.dat$method)
	 print(ggplot(dist.dat, aes(x = method,y=meanval,size=wts)) +
	 geom_point(color="blue")+
	 geom_segment(aes(x=c(1:dim(dist.dat)[1]),xend=c(1:dim(dist.dat)[1]),yend=meanval,y=0),size=1)+
	 geom_point(color="blue")+
	 scale_size_area(name  ="Weighting")+
	 scale_y_continuous(limits = c(0, NA))+
	 theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+
	 labs(x="Method",y="Natural Mortality"))
   }


 if(input$M_CV>0)
  {
  	dat.plot<-M.dists()
  	dat.plot$rowval<-factor(dat.plot$rowval,levels=unique(dat.plot$rowval))
  	
  	#col.dists<-colorRampPalette(c("#236192","#1D252D","#658D1B")) #Sounders colors
  	col.dists<-colorRampPalette(c("red", "yellow", "blue"))
  	#col.dists<-colorRampPalette(c("red", "yellow", "blue"))
  	print(ggplot(dat.plot, aes(x = dist,stat(count))) +
 		geom_density(aes(fill = factor(rowval)),alpha = 0.5)+
     	scale_fill_manual(values = col.dists(length(unique(dat.plot$rowval))),name="Method")+
 		#scale_fill_viridis(option="D",discrete=TRUE,name="Method")+
 		labs(x="Natural Mortality",y="Density"))
     	
  }

 	})

#Plot Composite M
 output$Mcomposite<- renderPlot({    
 

#    if(input$M_CV==0)
#   {
	pdf.Msamples<-M.dists()
	cdf.out<-ecdf(pdf.Msamples$dist)
  	Mcomposite.densityplot<-ggplot(data= pdf.Msamples,aes(dist))+
     	geom_density(fill="gray",bw="SJ")+
     	labs(x="Natural Mortality",y="Density")+ 
     	geom_vline(xintercept = quantile(cdf.out,0.5),color="darkblue",size=1.2)
 	print(Mcomposite.densityplot)
	
  #  #Calculate density function of point estimates	
  #  	M.densum<-density(M.sub.n0,weights=M.wts.sub.stand,cut=0)
  #  #Approximate the denisty function
  #  #f<- approxfun(M.densum$x, M.densum$y, yleft=0, yright=0)
  #  #Standardize densities
  #  	pdf_counts<-round(1000000*(M.densum$y/sum(M.densum$y)))
  #  #Expand densities to samples
  #  	pdf.samples<-unlist(mapply(rep,M.densum$x,pdf_counts))
  #  #Calculate the cdf
  #  	cdf.out<-ecdf(pdf.samples)
  #  #Plot the density function
  # 	M.densum.plot<- data.frame(x = M.densum$x, y = M.densum$y)
 	# Mcomposite.densityplot<- ggplot(data=M.densum.plot,aes(x,y,fill="blue"))+
  #    	geom_line(col="black")+
  #    	labs(x="Natural Mortality",y="Density")+ 
  #    	geom_area(fill="gray")+ 
  #    #scale_x_continuous(limits=c(0,quantile(M.densum$x,0.99999)))+
	 #   geom_vline(xintercept = quantile(cdf.out,0.5),color="darkblue",size=1.2)
	 #   print(Mcomposite.densityplot)
#   }
   
#if(input$M_CV>0)
#   {
#   	tot.samples<-1000000

#   }
   output$downloadMcompositedensityplot <- downloadHandler(
   filename = function() { paste0('Mcomposite_densityplot',Sys.time(), '.png')},
   content = function(file) {
     png(file, type='cairo',width=800,height=720)
     print(Mcomposite.densityplot)
     dev.off()},contentType = 'image/png') 
   output$downloadMcompositedist <- downloadHandler(
     filename = function() {  paste0("Mcomposite_samples",Sys.time(),".DMP") },
     content = function(file) {save(pdf.Msamples,file=file)}) 
    
   })
  }
#)

