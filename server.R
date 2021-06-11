require(shiny)
require(fishmethods)
require(ggplot2)
require(truncnorm)
require(data.table)
require(RColorBrewer) 
require(viridis)
require(reshape2)

#shinyServer(
  function(input, output) 
  {    
    
    Then_M<-function(Amax)
    {
      M_vals<-rep(NA,3)
      M_vals[1]<-4.889*Amax^-0.916 #Then_nls
      M_vals[2]<-exp(1.717-1.01*log(Amax)) #Then_lm
      #M_vals[2]<-5.109/Amax #Done in normal space, which is inappropriate
      M_vals[3]<-5.4/Amax #Hamel, which corrects the above Then by using the data in logspace
      return(M_vals)
    }
    
    M_ZM_AC<-function(k,Amax,t0)
    {
      Mvals.out<-c(NA,NA)
      Mvals.out[1]<- (3*k)/(exp(k*((0.302 * Amax)-t0))-1) 
      Mvals.out[2]<- (3*k)/(exp(k*((0.44 * Amax)-t0))-1)
      return(Mvals.out) 
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
    
    Hamel_M_k<-function(k)
    {
      M_val_Hamel_k<-k*1.753
      return(M_val_Hamel_k)
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

    Chen_N_Wat_Mage<-function(Age_in,k,t0)
    {
      if(anyNA(c(Age_in,k,t0))){M.out<-NA}
      else
      {
        a<-Age_in
        tM<--1/k*(log(abs(1-exp(k*t0)))+t0)
        a0<-1-exp(-k*(tM-t0))
        a1<-k*exp(-k*(tM-t0))
        a2<--0.5*k^2*exp(-k*(tM-t0))
        if(a<=tM){M.out<-k/(1-exp(-k*(a-t0)))}
        if(a>tM){M.out<-k/(a0+a1*(a-tM)+a2*(a-tM)^2)}
      }
      return(M.out)
    }

    Gislason_M_a<-function(Amax,Linf,k,t0)
    {
		Lts<-Linf*(1-exp(-k*(c(1:Amax)-t0)))
		Gis_Ms_a<-mapply(function(x) M.empirical(Linf=Linf,Kl=k,Bl=Lts[x],method=9)[1],x=1:length(Lts),SIMPLIFY=TRUE)
		return(Gis_Ms_a)
    }

  McCoyGillooly_M<-function(Mass,Temp)
  {
    McCGil_M<-((Mass/4)^-.25)*exp(-7540*((1/(273+Temp))-(1/293)))
    return(McCGil_M)
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
    
  #Calculate age-specific M values
   M_vals_ages<- reactive({
    CnW_M_a_VBGF<-Gislason_M_ages<-NA
    M_vals_ages<-data.table(Age=NA,CnW_M=CnW_M_a_VBGF,Gislason_M=Gislason_M_ages)
    #Calculate Gislason
    if(!(anyNA(c(input$Amax,input$Linf,input$k_vbgf,input$t0))))
    {Gislason_M_ages<-Gislason_M_a(input$Amax,input$Linf,input$k_vbgf,input$t0)}
    #Calculate Chen & Watanabe
    CnW_M_a_VBGF<-Chen_N_Wat_M(input$Amax,input$k_vbgf,input$t0,out.type = 0)
    #Put table together
    if(!is.na(input$Amax)){M_vals_ages<-data.table(Age=c(1:input$Amax),CnW_M=CnW_M_a_VBGF,Gislason_M=Gislason_M_ages)}
    #Create download object for the age-specific values csv file
    output$downloadCW_M_a <- downloadHandler(
     filename = function() {paste0("Age_specific_M_values", '.csv') },
     content = function(file) {write.csv(M_vals_ages, file=file)})  
   M_vals_ages
   })        
 
 #Calcualte individual estimates of M   
 M_vals_all<- reactive({
   fishlife.M.out<-Pauly80lt_M<-Pauly80wt_M<-AnC75_M<-Roff_M<-GnD_GSI_M<-PnW_M<-Lorenzen96_M<-Gislason_M<-McCGil_M<-NA
   #Fishlife
   if(length(strsplit(input$Genspp, " ")[[1]])>1)
   {
    fishlife.M.out<-try(fishlife.M(input$Genspp))
    if(inherits(fishlife.M.out,"try-error")==TRUE){fishlife.M.out<-NA}
  }
   #Longevity
   Then_M_Amax<-Then_M(input$Amax)
   #if(!(anyNA(c(input$k_vbgf,input$Amax)))){AnC75_M<-M.empirical(Kl=input$k_vbgf,tmax=input$Amax,method=4)[1]}
   ZMAC_M<-M_ZM_AC(k=input$k_vbgf,Amax=input$Amax,t0=input$t0)
   Chen_N_Wat_Ma<-Chen_N_Wat_Mage(input$Age_in,input$k_vbgf,input$t0)
   #VBGF
   Then_M_VBGF<-Then_VBGF(input$Linf,input$k_vbgf)
   Jensen_M_VBGF<-Jensen_M_k(input$k_vbgf) 
   Hamel_M_VBGF<-Hamel_M_k(input$k_vbgf) 
   if(!(anyNA(c(input$Linf,input$k_vbgf,input$Lt_in))))
   	{Gislason_M<-M.empirical(Linf=input$Linf,Kl=input$k_vbgf,Bl=input$Lt_in,method=9)[1]}
   #Maturity
   if(!(anyNA(c(input$k_vbgf,input$Amat)))){Roff_M<-M.empirical(Kl=input$k_vbgf,tm=input$Amat,method=5)[1]}
   Jensen_M_Amat<-Jensen_M_amat(input$Amat)
   Rikhter_Efanov_Amat<-Rikhter_Efanov_Amat_M(input$Amat)
   #Size-based
   if(!(anyNA(c(input$Wdry)))){PnW_M<-M.empirical(Wdry=input$Wdry,method=7)[1]}
   if(!(anyNA(c(input$Wwet)))){Lorenzen96_M<-M.empirical(Wwet=input$Wwet,method=8)[1]}
   if(!(anyNA(c(input$Linf,input$k_vbgf,input$Temp)))){Pauly80lt_M<-M.empirical(Linf=input$Linf,Kl=input$k_vbgf,TC=input$Temp,method=1)[1]}
   if(!(anyNA(c(input$Winf,input$kw,input$Temp)))){Pauly80wt_M<-M.empirical(Winf=input$Winf,Kw=input$kw,TC=input$Temp,method=2)[1]}
   if(!(anyNA(c(input$wdry,input$Temp)))){McCGil_M<-McCoyGillooly_M(input$Wdry,input$Temp)}
   #GSI
   #print(if(!(anyNA(c(input$GSI)))){GnD_GSI_M<-M.empirical(GSI=input$GSI,method=6)[1]})
   if(!(anyNA(c(input$GSI)))){GnD_GSI_M<-1.817*input$GSI} #Hamel update. Original value is 1.79
   #User inputs
   User_M<-as.numeric(trimws(unlist(strsplit(input$User_M,","))))
   if(length(User_M)==0)User_M<-NA
   M_users<-"User input"
   if(length(User_M)>1){M_users<-paste0("User input_",c(1:length(User_M)))}
   #Concatenate all M values
   M_vals_all<-c(fishlife.M.out,Then_M_Amax,Chen_N_Wat_Ma,ZMAC_M,Then_M_VBGF,Hamel_M_VBGF,Jensen_M_VBGF,Gislason_M,Pauly80lt_M,Roff_M,Jensen_M_Amat,Rikhter_Efanov_Amat,Pauly80wt_M,McCGil_M,PnW_M,Lorenzen96_M,GnD_GSI_M,User_M)
   M_methods<-c("FishLife","Then_nls","Then_lm","Hamel_Amax","Chen-Wat","ZM_CA_pel","ZM_CA_dem","Then_VBGF","Hamel_K","Jensen_K 1","Jensen_K 2","Gislason","Pauly_lt","Roff","Jensen_Amat","Ri_Ef_Amat","Pauly_wt","McC&Gil","PnW","Lorenzen","GSI",M_users)
   M_methods_vals_all<-data.table(Method=M_methods, M=M_vals_all)
   #Create object with all input parameter values
   M_parms_all<-c(input$M_CV,input$M_CV_type,input$Amax,input$Linf,input$k_vbgf,input$t0,input$Age_in,input$Lt_in,input$Amat,input$Temp,input$Winf,input$kw,input$Wdry,input$Wwet,input$GSI,User_M)
   M_parms_names<-c("CV","Error type","Max age","Linf","k","t0","Age","Length","Age mat","Temp","Winf","kW","Dry wt","Wet wt","GSI",M_users)
   M_parms_names_all<-data.table(Parameter=M_parms_names,Input=M_parms_all)
   M.out.parms.vals.ages<-list(Parameters= M_parms_names_all,M_estimates=M_methods_vals_all,Age_specific_values=M_vals_ages())
   #Create downlaod for M values, parameters and age-specific values
   output$downloadMandPs <- downloadHandler(
     filename = function() {  paste0("M_parms_values_byage_out",Sys.time(),".DMP") },
     content = function(file) {save(M.out.parms.vals.ages,file=file)}) 
   M_vals_all
   })


   output$Mplot <- renderPlot({
   M_vals_all<-M_vals_all()
   User_M<-as.numeric(trimws(unlist(strsplit(input$User_M,","))))
   M_users<-"User input"
   if(length(User_M)>1){M_users<-paste0("User input_",c(1:length(User_M)))}
   M_methods<-c("FishLife","Then_nls","Then_lm","Hamel_Amax","Chen-Wat","ZM_CA_pel","ZM_CA_dem","Then_VBGF","Hamel_k","Jensen_k 1","Jensen_k 2","Gislason","Pauly_lt","Roff","Jensen_Amat","Ri_Ef_Amat","Pauly_wt","McC&Gil","PnW","Lorenzen","GSI",M_users)
   M_types<-c("Meta-analysis",rep("Amax",3),rep("Amax:VBGF",3),rep("VBGF",5),rep("VBGF:Temp",1),"VBGF:Amat",rep("Amat",2),rep("Weight",4),rep("GSI",1),rep("User input",length(M_users)))
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
	  Mplots<-ggplot(M_vals_gg,aes(Method,as.numeric(as.character(M)),color=Input))+
           geom_point(size=4)+ylab("M")+xlab("Method")+
	  	     scale_y_continuous(limits = c(0, NA))+
	  	     #theme_minimal()+
           theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))

	  }
	 
      
	  if(input$M_CV>0 & input$M_CV_type=="lognormal")
	  {
      Mplots<-ggplot(M_vals_gg,aes(Method,as.numeric(as.character(M)),color=Input))+
           geom_point(size=4)+ylab("M")+xlab("Method")+
           scale_y_continuous(limits = c(0, NA))+
           theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+
   		   geom_pointrange(aes(ymin=qlnorm(0.025,log(as.numeric(as.character(M))),input$M_CV),ymax=qlnorm(0.975,log(as.numeric(as.character(M))),input$M_CV)))
      	   
      }

	if(input$M_CV>0 & input$M_CV_type=="normal")
	  {
      Mplots<-ggplot(M_vals_gg,aes(Method,as.numeric(as.character(M)),color=Input))+
           geom_point(size=4)+ylab("M")+xlab("Method")+
           theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+
   		     geom_pointrange(aes(ymin=qtruncnorm(0.025,a=0,mean=as.numeric(as.character(M)),sd=as.numeric(as.character(M))*input$M_CV),ymax=qtruncnorm(0.975,a=0,mean=as.numeric(as.character(M)),sd=as.numeric(as.character(M))*input$M_CV)))
      }
      print(Mplots)
   M_table<-data.frame(cbind(M_methods,M_vals_all))
   colnames(M_table)<-c("Method","M")
  # if(all(is.na(M_vals()))){return(NULL)}
   output$downloadMplots <- downloadHandler(
   filename = function() { paste0('Mplots',Sys.time(), '.png')},
   content = function(file) {
     png(file, type='cairo',width=800,height=720)
     print(Mplots)
     dev.off()},contentType = 'image/png') 
   output$downloadMs <- downloadHandler(
    filename = function() {paste0("M_values", '.csv') },
    content = function(file) {write.csv(M_table, file=file)}
  )
 })

 output$Mplot_ages <- renderPlot({
    
    M_vals_ages_melt<-data.table(melt(M_vals_ages(),id.vars=1,variable.name="Method",value.name="M"))
    Mplot_ages<-ggplot(M_vals_ages_melt,aes(Age,M,color=Method))+
    geom_point(size=4)+
    #scale_y_continuous(limits = c(0, NA))+
    labs(x="Age",y="Natural mortality")
    print(Mplot_ages)
     output$downloadMplot_ages <- downloadHandler(
    filename = function() { paste0('Mplot_ages',Sys.time(), '.png')},
    content = function(file) {
     png(file, type='cairo',width=800,height=720)
     print(Mplot_ages)
     dev.off()},contentType = 'image/png') 
  
 })

# Show the first "n" observations
 output$Mtable <- renderTable({
   fishlife.M.out<-AnC75_M<-Gislason_M<-NA
   if(length(strsplit(input$Genspp, " ")[[1]])>1)
   {
    fishlife.M.out<-try(fishlife.M(input$Genspp))
    if(inherits(fishlife.M.out,"try-error")==TRUE){fishlife.M.out<-NA}
  }
   Then_M_Amax<-Then_M(input$Amax)
   #if(!(anyNA(c(input$k_vbgf,input$Amax)))){AnC75_M<-M.empirical(Kl=input$k_vbgf,tmax=input$Amax,method=4)[1]}
   #if(!(anyNA(c(input$k_vbgf,input$Amax,input$t0)))){ZMAC_M<-M_ZM_AC(k=input$k_vbgf,Amax=input$Amax,t0=input$t0)}
   ZMAC_M<-M_ZM_AC(k=input$k_vbgf,Amax=input$Amax,t0=input$t0)
   Then_M_VBGF<-Then_VBGF(input$Linf,input$k_vbgf)
   Jensen_M_VBGF<-Jensen_M_k(input$k_vbgf) 
   Hamel_M_VBGF<-Hamel_M_k(input$k_vbgf) 
   if(!(anyNA(c(input$Linf,input$k_vbgf,input$Lt_in)))){Gislason_M<-M.empirical(Linf=input$Linf,Kl=input$k_vbgf,Bl=input$Lt_in,method=9)[1]}
   Chen_N_Wat_Ma<-Chen_N_Wat_Mage(input$Age_in,input$k_vbgf,input$t0)
   
   M_vals_all<-c(fishlife.M.out,Then_M_Amax,Chen_N_Wat_Ma,ZMAC_M,Then_M_VBGF,Hamel_M_VBGF,Jensen_M_VBGF,Gislason_M)
   M_methods<-c("FishLife","Then_nls","Then_lm","Hamel_Amax","Chen-Wat","ZM_CA_pel","ZM_CA_dem","Then_VBGF","Hamel_k","Jensen_k 1","Jensen_k 2","Gislason")
   M_table<-data.frame(cbind(M_methods,signif(M_vals_all,3)))
   colnames(M_table)<-c("Method","M")
   #rownames(M_table)<-M_methods
   M_table
  })
# Show the first "n" observations
 output$Mtable2 <- renderTable({
   Jensen_M_Amat<-Pauly80lt_M<-Pauly80wt_M<-Roff_M<-GnD_GSI_M<-PnW_M<-Lorenzen96_M<-Gislason_M<-McCGil_M<-NA
   if(!(anyNA(c(input$k_vbgf,input$Amat)))){Roff_M<-M.empirical(Kl=input$k_vbgf,tm=input$Amat,method=5)[1]}
   Jensen_M_Amat<-Jensen_M_amat(input$Amat)
   Rikhter_Efanov_Amat<-Rikhter_Efanov_Amat_M(input$Amat)
   if(!(anyNA(c(input$Wdry)))){PnW_M<-M.empirical(Wdry=input$Wdry,method=7)[1]}
   if(!(anyNA(c(input$Wwet)))){Lorenzen96_M<-M.empirical(Wwet=input$Wwet,method=8)[1]}
   if(!(anyNA(c(input$Linf,input$k_vbgf,input$Temp)))){Pauly80lt_M<-M.empirical(Linf=input$Linf,Kl=input$k_vbgf,TC=input$Temp,method=1)[1]}
   if(!(anyNA(c(input$Winf,input$kw,input$Temp)))){Pauly80wt_M<-M.empirical(Winf=input$Winf,Kw=input$kw,TC=input$Temp,method=2)[1]}
   if(!(anyNA(c(input$GSI)))){GnD_GSI_M<-1.79*input$GSI}
   if(!(anyNA(c(input$Wdry,input$Temp)))){McCGil_M<-McCoyGillooly_M(input$Wdry,input$Temp)}
   #if(!(anyNA(c(input$GSI)))){GnD_GSI_M<-M.empirical(GSI=input$GSI,method=6)[1]}
   #User_M<-input$User_M
   
   M_vals_all<-c(Pauly80lt_M,Roff_M,Jensen_M_Amat,Rikhter_Efanov_Amat,Pauly80wt_M,McCGil_M,PnW_M,Lorenzen96_M,GnD_GSI_M)
   M_methods<-c("Pauly_lt","Roff","Jensen_Amat","Ri_Ef_Amat","Pauly_wt","McC&Gil","PnW","Lorenzen","GSI")
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
   			input$Then_nls,
   			input$Then_lm,
   			input$Hamel_Amax,
   			input$Chen_Wat,
   			input$ZM_CA_pel,
   			input$ZM_CA_dem,
        input$Then_VBGF,
   			input$Hamel_VBGF,
        input$Jensen_VBGF_1,
   			input$Jensen_VBGF_2,
   			input$Gislason,
   			input$Pauly_lt,
   			input$Roff,
   			input$Jensen_Amat,
   			input$Ri_Ef_Amat,
   			input$Pauly_wt,
   			input$McGl,
   			input$PnW,
   			input$Lorenzen,
   			input$Gonosoma,
   			rep(input$UserM_wt,length(User_M)))
   names(M.wts)<-c("FishLife",
   					"Then_nls",
   					"Then_lm",
   					"Hamel_Amax",
   					"Chen-Wat",
   					"ZM_CA_pel",
   					"ZM_CA_dem",
            "Then_VBGF",
   					"Hamel_k",
            "Jensen_k1",
   					"Jensen_k2",
   					"Gislason",
   					"Pauly_lt",
   					"Roff",
   					"Jensen_Amat",
   					"Ri_Ef_Amat",
   					"Pauly_wt",
   					"McC&Gil",
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
   samp.num<-input$samp.num
   samps<-samp.num*M.wts.sub.stand
	
	M.CV.method <- data.table(meanM = M.sub.n0,
                       CVval = input$M_CV,
                       Samples = samps,
                       Original.wts = M.wts.sub.n0,
                       Stand.Weights=M.wts.sub.stand,
                       Method= names(M.wts.sub.stand))
	M.CV.method
	}
	})

 	M.dists<- reactive({
		if(input$M_CV_type=="lognormal")
		{
			M.dists <- rbindlist(lapply(1:dim(M.CV.method())[1], 
                        function(x) data.table(Method = M.CV.method()$Method[x], 
                                               Mval = rlnorm(M.CV.method()[x, Samples],
                                               log(M.CV.method()[x, meanM]), 
                                               M.CV.method()[x, CVval]))))
		}
		if(input$M_CV_type=="normal")
		{
			M.dists <- rbindlist(lapply(1:dim(M.CV.method())[1], 
                        function(x) data.table(Method = M.CV.method()$Method[x], 
                                               Mval = rtruncnorm(M.CV.method()[x, Samples],
                                               a=0,
                                               b=Inf, 
                                               M.CV.method()[x, meanM], 
                                               M.CV.method()[x, meanM*CVval]))))
		}
  M.dists
 	})

#Plot Individual distributions
 output$Mdistplots<- renderPlot({ 
   req(!all(is.na(M_vals_all())))
   dist.dat<-M.CV.method()
   if(input$M_CV==0)
   {
	 dist.dat$Method<-factor(dist.dat$Method,levels=dist.dat$Method)
	 Mdist_plots<-ggplot(dist.dat, aes(x = Method,y=meanM,size=Stand.Weights)) +
	 geom_point(color="blue")+
	 geom_segment(aes(x=c(1:dim(dist.dat)[1]),xend=c(1:dim(dist.dat)[1]),yend=meanM,y=0),size=1)+
	 geom_point(color="blue")+
	 scale_size_area(name  ="Weighting")+
	 scale_y_continuous(limits = c(0, NA))+
	 theme(axis.text.x = element_text(angle = 90, hjust = 1,vjust=0.5))+
	 labs(x="Method",y="Natural Mortality")
   print(Mdist_plots)
   }


 if(input$M_CV>0)
  {
  	dat.plot<-M.dists()
  	dat.plot$Method<-factor(dat.plot$Method,levels=unique(dat.plot$Method))
  	#col.dists<-colorRampPalette(c("#236192","#1D252D","#658D1B")) #Sounders colors
  	col.dists<-colorRampPalette(c("red", "yellow", "blue"))
  	#col.dists<-colorRampPalette(c("red", "yellow", "blue"))
  	Mdist_plots<-ggplot(dat.plot, aes(x = Mval,fill = factor(Method))) +
 		#geom_density(aes(fill = factor(Method),alpha = 0.5))+
    geom_area(stat="bin",binwidth=0.005,alpha = 0.5)+
     	scale_fill_manual(values = col.dists(length(unique(dat.plot$Method))),name="Method")+
 		#scale_fill_viridis(option="D",discrete=TRUE,name="Method")+
 		theme_minimal()+
    labs(x="Natural Mortality",y="Density")
    print(Mdist_plots)
  }
   Mdist_values<-dist.dat
   output$downloadMdensityplots <- downloadHandler(
   filename = function() { paste0('Mdensityplots',Sys.time(), '.png')},
   content = function(file) {
     png(file, type='cairo',width=800,height=720)
     print(Mdist_plots)
     dev.off()},contentType = 'image/png') 
   output$downloadMdistvals <- downloadHandler(
     filename = function() {  paste0("Mdist_values",Sys.time(),".DMP") },
     content = function(file) {save(Mdist_values,file=file)}) 
 	})

#Calculate updated M prior given changes in bandwidth
priorMupdate<-reactive({
    Msamples<-M.dists()
    Mvals<-Msamples$Mval
    Mvals.update<-density(Mvals,adjust=input$ad.bw,from=0)  
    priorMupdate <- approx(
    cumsum(Mvals.update$y)/sum(Mvals.update$y),
    Mvals.update$x,
    runif(input$samp.num)
    )$y
    priorMupdate
  })

#Plot Composite M
 output$Mcomposite<- renderPlot({    
 	req(!all(is.na(M_vals_all())))
  Msamples<-M.dists()
  priorMupdate<-priorMupdate()
  cdf.out<-ecdf(Msamples$Mval)
  	Mcomposite.densityplot<-ggplot(data= Msamples,aes(Mval))+
     	geom_density(fill="gray",bw="SJ",adjust=input$ad.bw)+
     	labs(x="Natural Mortality",y="Density")+ 
     	geom_vline(xintercept = quantile(cdf.out,0.5),color="darkblue",size=1.2)+
      theme_minimal()
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
   
   output$downloadMcompositedensityplot <- downloadHandler(
   filename = function() { paste0('Mcomposite_densityplot',Sys.time(), '.png')},
   content = function(file) {
     png(file, type='cairo',width=800,height=720)
     print(Mcomposite.densityplot)
     dev.off()},contentType = 'image/png') 
   output$downloadMcompositedist <- downloadHandler(
     filename = function() {  paste0("Mcomposite_samples",Sys.time(),".DMP") },
     content = function(file) {save(Msamples,file=file)})     
    output$downloadMcompositedistupdated <- downloadHandler(
     filename = function() { paste0("UpdatedMcomposite_samples",Sys.time(),".DMP") },
     content = function(file) {save(priorMupdate,file=file)})  
   })

  }
#)

