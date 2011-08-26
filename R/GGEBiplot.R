GGEBiplot <-
function ()
{
    library(RODBC)
    require(tcltk)
    library(tkrplot)
    library(rgl)
    tclRequire("BWidget")
    #
    # Variables
    #
    optioncentering <- "2.Tester-Centered G+GE"
    optionscaling  <- "0.No scaling"
    optionSVP <- "GH -(Column Metric Preserving)"
    datascaling <- c("0.No scaling","1.Std Deviation (SD)")            
    datacentering <- c("0.No centering","1.Global-Centered E+G+GE","2.Tester-Centered G+GE","3.Double-Centered GE")    
    dataSVP <- c("JK -(Row Metric Preserving)","GH -(Column Metric Preserving)","HJ -(Dual Metric Preserving)","SQ - Symmetrical")            
    wintitle <- "GGE Biplot"    
    coltitle <- "black"
    background <- "white"    
    colenv <- "blue"
    colgenotype <- "green"
    centro <- c(0, 0)    
    symbol = NA_integer_
    subtitle <- NULL
    ejes <- array()
    showtitle <- tclVar("1")
    showboth <- tclVar("0") 
    vaxis <- tclVar("0")
    showguidelines <- tclVar("1")
    showcircles <- tclVar("0")
    scaling <- tclVar("0")
    centering <- tclVar("2")
    svp <- tclVar("1")
    vrank <- tclVar("1")
    TypeGraph <- 1 #Biplot
    matrixdata <- NULL # Matriz de datos
    coordgenotype <- NULL # Coordenadas de los genotipos
    coordenviroment <- NULL # Coordenadas de los ambientes
    labelenv <- NULL # Etiquetas de los ambientes
    labelgen <- NULL # Etiquetas de los genotipos
    venvironment <- -1
    vgenotype <- -1
    vgenotype1 <- -1
    vgenotype2 <- -1
    dimension1 <- 1
    dimension2 <- 2
    #
    # Función que abre una ventana de dialogo.
    #
    modalDialog <- function(title,question,entryInit,entryWidth=20,returnValOnCancel="ID_CANCEL")
      {
        dlg <- tktoplevel()
        tkwm.deiconify(dlg)
        tkgrab.set(dlg)
        tkfocus(dlg)
        tkwm.title(dlg,title)
        textEntryVarTcl <- tclVar(paste(entryInit))
        textEntryWidget <- tkentry(dlg,width=paste(entryWidth),textvariable=textEntryVarTcl)
        tkgrid(tklabel(dlg,text="       "))
        tkgrid(tklabel(dlg,text=question),textEntryWidget)
        tkgrid(tklabel(dlg,text="       "))
        ReturnVal <- returnValOnCancel
        onOK <- function()
        {
          ReturnVal <<- tclvalue(textEntryVarTcl)
          tkgrab.release(dlg)
          tkdestroy(dlg)
        }
        onCancel <- function()
        {
          ReturnVal <<- returnValOnCancel
          tkgrab.release(dlg)
          tkdestroy(dlg)
        }
        OK.but     <-tkbutton(dlg,text="   OK   ",command=onOK)
        Cancel.but <-tkbutton(dlg,text=" Cancel ",command=onCancel)
        tkgrid(OK.but,Cancel.but)
        tkgrid(tklabel(dlg,text="    "))
        tkfocus(dlg)
        tkwait.window(dlg)
        return(ReturnVal)
    }
    #
    # Función cambia color
    #
    ChangeColorv <- function(color)
    {
      colorv <<- tclvalue(tcl("tk_chooseColor", initialcolor = color,
      title = "Choose a color"))
      if (nchar(colorv) > 0) return(colorv)
    }  
    #
    # Función modelos
    # 
    Models <- function ()        
    {
      # Carga la matriz de datos, los nombres de las filas y las columnas.
      labelenv <<- array(,ncol(mdata)-1)
      labelgen <<- array(,nrow(mdata))        
      for (i in 1:ncol(mdata)-1) labelenv[i] <<- colnames(mdata)[i+1]
      labelgen <<- mdata[,1]
      matrixdata <<- matrix(,nrow(mdata),ncol(mdata)-1)        
      for (i in 1:nrow(matrixdata))         
        for (j in 1:ncol(matrixdata)) 
         matrixdata[i,j] <<- mdata[i,j+1]        
      for (i in 1:ncol(diag(svd(matrixdata)$d))) 
        ejes[i] <<- paste("AXIS", i, sep = "")                                                 
      # Opción de centrado             
      if (optioncentering == "0.No centering")
        {
          centering <<- tclVar("0")
        }
      if (optioncentering == "1.Global-Centered E+G+GE")
        {
          meanData <<- mean(matrixdata)
          matrixdata <<- matrixdata-meanData
          centering <<- tclVar("1")
        }
      if (optioncentering == "2.Tester-Centered G+GE")
        {
          meancolData <<- colMeans(matrixdata)
          for (i in  1:nrow(matrixdata))
            for (j in 1:ncol(matrixdata))
              matrixdata[i,j] <<- matrixdata[i,j]-meancolData[j]
          centering <<- tclVar("2")    
        }
      if (optioncentering == "3.Double-Centered GE")
        {
          meanData <<- mean(matrixdata)
          meancolData <<- colMeans(matrixdata)
          meanrowData <<- rowMeans (matrixdata)
          for (i in  1:nrow(matrixdata))
            for (j in 1:ncol(matrixdata))
              matrixdata[i,j] <<- matrixdata[i,j]+ meanData - meancolData[j] - meanrowData[i]
          centering <<- tclVar("3")
        }
      # Opción de escalado
      if (optionscaling == "0.No scaling")
        {
          scaling <<- tclVar("0")
        }
      if (optionscaling == "1.Std Deviation (SD)")
        {
          scaling <<- tclVar("1")
          desviation <<- array(,dim=ncol(matrixdata))
          for (j in 1:ncol(matrixdata)) desviation[j] <<- sqrt(var(matrixdata[,j]))
          for (i in  1:nrow(matrixdata)) 
            for (j in 1:ncol(matrixdata))
              matrixdata[i,j] <<- matrixdata[i,j]/desviation[j]
        }
      # Opción de SVP
      if (optionSVP == "JK -(Row Metric Preserving)")
        {
          coordgenotype <<- svd(matrixdata)$u %*% diag(svd(matrixdata)$d)
          coordenviroment <<- svd(matrixdata)$v
          d1 <<- (max(coordenviroment[,dimension1])-min(coordenviroment[,dimension1]))/(max(coordgenotype[,dimension1])-min(coordgenotype[,dimension1]))
          d2 <<- (max(coordenviroment[,dimension2])-min(coordenviroment[,dimension2]))/(max(coordgenotype[,dimension2])-min(coordgenotype[,dimension2]))
          d <<- max(d1,d2)
          coordenviroment <<- coordenviroment/d
          svp <<- tclVar("0")
        }
      if (optionSVP == "GH -(Column Metric Preserving)")
        {
          coordgenotype <<- svd(matrixdata)$u
          coordenviroment <<- svd(matrixdata)$v %*% diag(svd(matrixdata)$d)
          d1 <<- (max(coordgenotype[,dimension1])-min(coordgenotype[,dimension1]))/(max(coordenviroment[,dimension1])-min(coordenviroment[,dimension1]))
          d2 <<- (max(coordgenotype[,dimension2])-min(coordgenotype[,dimension2]))/(max(coordenviroment[,dimension2])-min(coordenviroment[,dimension2]))
          d <<-  max(d1,d2)
          coordgenotype <<- coordgenotype/d              
          svp <<- tclVar("1")
        }
      if (optionSVP == "SQ - Symmetrical")
        {
          coordgenotype <<- svd(matrixdata)$u %*% diag(sqrt(svd(matrixdata)$d))
          coordenviroment <<- svd(matrixdata)$v%*% diag(sqrt(svd(matrixdata)$d))
          svp <<- tclVar("3")
        }
      if (optionSVP == "HJ -(Dual Metric Preserving)")
        {
          coordgenotype <<- svd(matrixdata)$u %*% diag(svd(matrixdata)$d)
          coordenviroment <<- svd(matrixdata)$v%*% diag(svd(matrixdata)$d)                        
          svp <<- tclVar("2")
        }
    }
    #
    # Función carga fichero
    #  
    Addfile <- function()
      {                        
        valorespropios <<- svd(matrixdata)$d
        vartotal <<- round(as.numeric(sum(valorespropios^2)),2)
        varexpl <<- round(as.numeric((valorespropios^2/vartotal) *100),2)                                                   
        genfile <<- as.data.frame(coordgenotype[,dimension1:dimension2])
        rownames(genfile) <<- labelgen
        colnames(genfile) <<- ejes[dimension1:dimension2]                   
        envfile <<- as.data.frame(coordenviroment[,dimension1:dimension2])              
        rownames(envfile) <<- labelenv
        colnames(envfile) <<- ejes[dimension1:dimension2]                                
        coordgencuad <<- coordgenotype^2
        CRFqEi <<- coordgencuad
        sumacuagen <- rowSums(coordgencuad)
        CRFqEi[, 1] <<- round(((coordgencuad)[, dimension1] * 1000)/sumacuagen,0)
        CRFqEi[, 2] <<- round(((coordgencuad)[, dimension2] * 1000)/sumacuagen,0)
        CRFqEi <<- as.data.frame(CRFqEi[,dimension1:dimension2])
        rownames(CRFqEi) <<- labelgen
        colnames(CRFqEi) <<- ejes[dimension1:dimension2]
        coordenvcuad <<- coordenviroment^2              
        CRFqEj <<- coordenvcuad
        sumacuaenv <<- rowSums(coordenvcuad)
        CRFqEj[, 1] <<- round(((coordenvcuad)[, dimension1] * 1000)/(sumacuaenv),0)
        CRFqEj[, 2] <<- round(((coordenvcuad)[, dimension2] * 1000)/(sumacuaenv),0)
        CRFqEj <<- as.data.frame(CRFqEj[,1:2])
        rownames(CRFqEj) <<- labelenv
        colnames(CRFqEj) <<- ejes[dimension1:dimension2]                              
        cat("GGE BIPLOT",file = "Results1.xls")
        file.append("Results1.xls","temp.xls")
        cat("\n",file="temp.xls")
        file.append("Results1.xls","temp.xls")
        cat("Centered by: ",file="temp.xls")
        file.append("Results1.xls","temp.xls")
        cat(optioncentering,file="temp.xls")
        file.append("Results1.xls","temp.xls")
        cat("\n",file="temp.xls")             
        file.append("Results1.xls","temp.xls")
        cat("Scaled (Divided) by: ",file="temp.xls")
        file.append("Results1.xls","temp.xls")
        cat(optionscaling,file="temp.xls")
        file.append("Results1.xls","temp.xls")
        cat("\n",file="temp.xls")     
        file.append("Results1.xls","temp.xls")                              
        cat("SVP: ",file="temp.xls")
        file.append("Results1.xls","temp.xls")
        cat(optionSVP,file="temp.xls")
        file.append("Results1.xls","temp.xls") 
        cat("\n",file="temp.xls")             
        file.append("Results1.xls","temp.xls")              
        cat("\n",file="temp.xls")             
        file.append("Results1.xls","temp.xls")                                     
        cat("Eigenvalues and variance explained",file="temp.xls")       
        file.append("Results1.xls","temp.xls")                   
        write.table(round(svd(matrixdata)$d,3),file="temp.xls",sep="\t",dec=",")        
        file.append("Results1.xls","temp.xls")                                         
        cat("\n",file="temp.xls")
        file.append("Results1.xls","temp.xls")                                 
        cat("Row coordinates:",file="temp.xls")       
        file.append("Results1.xls","temp.xls")
        cat("\n",file="temp.xls")
        file.append("Results1.xls","temp.xls")                                                                               
        write.table(round(genfile,3),file="temp.xls",sep="\t",dec=",")        
        file.append("Results1.xls","temp.xls")                                         
        cat("\n",file="temp.xls")
        file.append("Results1.xls","temp.xls")                                 
        cat("Column coordinates:",file="temp.xls")       
        file.append("Results1.xls","temp.xls")     
        cat("\n",file="temp.xls")
        file.append("Results1.xls","temp.xls")                                                                                        
        write.table(round(envfile,3),file="temp.xls",sep="\t",dec=",")        
        file.append("Results1.xls","temp.xls")  
        cat("\n",file="temp.xls")
        file.append("Results1.xls","temp.xls")                                                                                        
        cat("\n",file="temp.xls")
        file.append("Results1.xls","temp.xls")                                                                                                                    
        cat("RELATIVE CONTRIBUTIONS OF THE FACTOR TO THE ELEMENT:",file="temp.xls")       
        file.append("Results1.xls","temp.xls") 
        cat("\n",file="temp.xls")
        file.append("Results1.xls","temp.xls")                                                                                                                                  
        cat("Row Contributions ----------",file="temp.xls")       
        file.append("Results1.xls","temp.xls")                                                                                                                                                                                                                                                                              
        cat("\n",file="temp.xls")
        file.append("Results1.xls","temp.xls")                                                                                                                                                                                                                                                                                            
        write.table(CRFqEi,file="temp.xls",sep="\t",dec=",")        
        file.append("Results1.xls","temp.xls")                
        cat("\n",file="temp.xls")
        file.append("Results1.xls","temp.xls")                                                                                                                                                                                                                                                                                                          
        cat("Column Contributions ----------",file="temp.xls")       
        file.append("Results1.xls","temp.xls")                                                                                                                                                                                                                                                                              
        cat("\n",file="temp.xls")
        file.append("Results1.xls","temp.xls")                                                                                                                                                                                                                                                                                            
        write.table(CRFqEj,file="temp.xls",sep="\t",dec=",")        
        file.append("Results1.xls","temp.xls")                              
        cat("\n", file = "temp.xls")                  
        file.append("Results1.xls","temp.xls")        
        file.show("Results1.xls")        
        file.remove("temp.xls")
    }                          
    #
    # Función plot
    #
    plotFunctiond <- function(screen = TRUE) 
    {
      valorespropios <<- svd(matrixdata)$d
      vartotal <<- round(as.numeric(sum(valorespropios^2)),2)
      varexpl <<- round(as.numeric((valorespropios^2/vartotal) *100),2)                                                       
      params <- par(bg = background)
      plot(rbind(coordgenotype,coordenviroment), main = wintitle, type = "n", asp=1,col.main = coltitle, xlab=paste(ejes[dimension1],varexpl[dimension1],"%",sep=" ",sub=subtitle)
           ,ylab=paste(ejes[dimension2],varexpl[dimension2],"%",sep=" "))
      if (tclvalue(showguidelines)== "1") abline (h=0,v=0,lty="dotted")                   
      if (TypeGraph == 1) # Biplot
        {
          if (tclvalue(showboth)== "0" || tclvalue(showboth)== "1")
            { 
              points(coordgenotype[,dimension1], coordgenotype[,dimension2],pch = symbol, col = colgenotype)
              text(coordgenotype[,dimension1],coordgenotype[,dimension2],labels=labelgen,col=colgenotype,cex=.8)
            }
          if (tclvalue(showboth)== "0" || tclvalue(showboth)== "2")
            {
              arrows(centro[1], centro[2], coordenviroment[,dimension1], coordenviroment[,dimension2], col = colenv,lty = "dotted", length = 0.05)
              points(centro[1], centro[2], pch = 18, col = "black")              
              text(coordenviroment[,dimension1], coordenviroment[,dimension2],labels=labelenv,col=colenv,cex=1)
            }
        }  
      if (TypeGraph == 2) # Evalua enviroment
      {
          points(coordgenotype[,dimension1], coordgenotype[,dimension2],pch = symbol, col = colgenotype)
          text(coordgenotype[,dimension1], coordgenotype[,dimension2],labels=labelgen,col=colgenotype,cex=.8)
          abline(a=0,b=coordenviroment[venvironment,dimension2]/coordenviroment[venvironment,dimension1],col="red",lty="solid")
          abline(a=0,b=-coordenviroment[venvironment,dimension1]/coordenviroment[venvironment,dimension2],col="red",lty="solid")                
          arrows(centro[1], centro[2], coordenviroment[venvironment,dimension1], coordenviroment[venvironment,dimension2], col = "red", lty = "solid", length = 0.1)
          text(coordenviroment[venvironment,dimension1],coordenviroment[venvironment,dimension2],labels=labelenv[venvironment],col="red")                
          for (i in 1:nrow(matrixdata))
          {
              x <- solve(matrix(c(-coordenviroment[venvironment,dimension2],coordenviroment[venvironment,dimension1],coordenviroment[venvironment,dimension1],coordenviroment[venvironment,dimension2]),nrow=2)
                   ,matrix(c(0,coordenviroment[venvironment,dimension1]*coordgenotype[i,dimension1]+coordenviroment[venvironment,dimension2]*coordgenotype[i,dimension2]),ncol=1))
              segments(coordgenotype[i,dimension1],coordgenotype[i,dimension2],x[1],x[2],lty="dotted")
          }
      }
      if (TypeGraph == 3) # Evalua genotypes
      {
         points(coordenviroment[,dimension1], coordenviroment[,dimension2],pch = symbol, col = colenv)
         text(coordenviroment[,dimension1], coordenviroment[,dimension2],labels=labelenv,col=colenv,cex=.8)
         abline(a=0,b=coordgenotype[vgenotype,dimension2]/coordgenotype[vgenotype,dimension1],col="red",lty="solid")
         abline(a=0,b=-coordgenotype[vgenotype,dimension1]/coordgenotype[vgenotype,dimension2],col="red",lty="solid")                
         arrows(centro[1], centro[2], coordgenotype[vgenotype,dimension1], coordgenotype[vgenotype,dimension2], col = "red", lty = "solid", length = 0.1)
         text(coordgenotype[vgenotype,dimension1],coordgenotype[vgenotype,dimension2],labels=labelgen[vgenotype],col="red")                
         for (i in 1:ncol(matrixdata))
           {
             x <- solve(matrix(c(-coordgenotype[vgenotype,dimension2],coordgenotype[vgenotype,dimension1],coordgenotype[vgenotype,dimension1],coordgenotype[vgenotype,dimension2]),nrow=2)
                 ,matrix(c(0,coordgenotype[vgenotype,dimension1]*coordenviroment[i,dimension1]+coordgenotype[vgenotype,dimension2]*coordenviroment[i,dimension2]),ncol=1))
             segments(coordenviroment[i,dimension1],coordenviroment[i,dimension2],x[1],x[2],lty="dotted")
           }                
      }
      if (TypeGraph == 4) # Relation among environments
      {
        text(coordgenotype[,dimension1],coordgenotype[,dimension2],labels=labelgen,col=colgenotype,cex=.8)            
        arrows(centro[1], centro[2], coordenviroment[,dimension1], coordenviroment[,dimension2], col = colenv,lty = "solid", length = 0.05)
        points(centro[1], centro[2], pch = 18, col = "black")              
        text(coordenviroment[,dimension1],coordenviroment[,dimension2],labels=labelenv,col=colenv,cex=1)
        if (tclvalue(showcircles) == "1")
        {
          radio <<- max((max(coordenviroment[dimension1,])-min(coordenviroment[dimension1,])),(max(coordenviroment[dimension2,])-min(coordenviroment[dimension2,])))/10
          for (i in 1:5) symbols(0,0,circles=radio*i,add = TRUE, inches = FALSE,fg = "black")          
        }       
      }
      if (TypeGraph == 5) # Compare genotypes
      {
         text(coordenviroment[,dimension1],coordenviroment[,dimension2],labels=labelenv,col=colenv)                         
         text(coordgenotype[,dimension1],coordgenotype[,dimension2],labels=labelgen,col=colgenotype)                
         symbols(coordgenotype[vgenotype1,dimension1],coordgenotype[vgenotype1,dimension2],circles=0.2,add = TRUE, inches = FALSE,fg = colgenotype)                           
         symbols(coordgenotype[vgenotype2,dimension1],coordgenotype[vgenotype2,dimension2],circles=0.2,add = TRUE, inches = FALSE,fg = colgenotype)                                    
         segments(coordgenotype[vgenotype1,dimension1],coordgenotype[vgenotype1,dimension2],
                  coordgenotype[vgenotype2,dimension1],coordgenotype[vgenotype2,dimension2],col="red",lty="solid")
         abline(a=0,b=-(coordgenotype[vgenotype1,dimension1]-coordgenotype[vgenotype2,dimension1])
                /(coordgenotype[vgenotype1,dimension2]- coordgenotype[vgenotype2,dimension2]),col="red",lty="solid")                         
      } 
      if (TypeGraph == 6) # Polygon view
      {
        points(coordgenotype[,dimension1], coordgenotype[,dimension2],pch = symbol, col = colgenotype)
        text(coordgenotype[,dimension1], coordgenotype[,dimension2],labels=labelgen,col=colgenotype,cex=0.8)
        points(centro[1], centro[2], pch = 18, col = "black")              
        text(coordenviroment[,dimension1],coordenviroment[,dimension2],labels=labelenv,col=colenv,cex=0.8)
        indice <<- c(chull(coordgenotype[,dimension1],coordgenotype[,dimension2]))
        polygon (coordgenotype[indice,dimension1],coordgenotype[indice,dimension2],border="black")
        i <<- 1
        while (is.na(indice[i+1])== FALSE)
        {
          abline(a=0,b=-(coordgenotype[indice[i],dimension1]-coordgenotype[indice[i+1],dimension1])
               /(coordgenotype[indice[i],dimension2]- coordgenotype[indice[i+1],dimension2]),col="red",lty="solid")                                 
          i <<- i+1               
        }
        abline(a=0,b=-(coordgenotype[indice[i],dimension1]-coordgenotype[indice[1],dimension1])
               /(coordgenotype[indice[i],dimension2]- coordgenotype[indice[1],dimension2]),col="red",lty="solid")                                         
      }       
      if (TypeGraph == 7) # Discrimitiveness vs. representativenss     
      {
        text(coordgenotype[,dimension1],coordgenotype[,dimension2],labels=labelgen,col=colgenotype,cex=.8)      
        segments(centro[1], centro[2], coordenviroment[,dimension1], coordenviroment[,dimension2], col = colenv,lty = "dotted")
        points(centro[1], centro[2], pch = 18, col = "black")              
        text(coordenviroment[,dimension1],coordenviroment[,dimension2],labels=labelenv,col=colenv,cex=1)
        arrows(centro[1], centro[2], mean(coordenviroment[,dimension1]), mean(coordenviroment[,dimension2]), col = colenv,lty = "solid", length = 0.1)
        symbols(mean(coordenviroment[,dimension1]),mean(coordenviroment[,dimension2]),circles=0.1,add = TRUE, inches = FALSE,fg = colenv)                  
        abline(a=0,b=mean(coordenviroment[,dimension2])/mean(coordenviroment[,dimension1]),col=colenv,lty="solid",cex=2)                                                 
        radio <<- max((max(coordenviroment[dimension1,])-min(coordenviroment[dimension1,])),(max(coordenviroment[dimension2,])-min(coordenviroment[dimension2,])))/10
        for (i in 1:5) symbols(0,0,circles=radio*i,add = TRUE, inches = FALSE,fg = "black")                  
      }
      if (TypeGraph == 9) # Mean vs. stability
      {
        med1 <<- mean(coordenviroment[,dimension1])
        med2 <<- mean(coordenviroment[,dimension2])        
        abline(a=0,b=med2/med1,col=colgenotype,lty="solid",cex=2) 
        abline(a=0,b=-med1/med2,col=colgenotype,lty="solid",cex=2)         
        arrows(centro[1], centro[2], med1, med2, col = colgenotype,lty = "solid", length = 0.1)                                                                
        text(coordgenotype[,dimension1],coordgenotype[,dimension2],labels=labelgen,col=colgenotype,cex=.8)      
        text(coordenviroment[,dimension1],coordenviroment[,dimension2],labels=labelenv,col=colenv,cex=1)        
        symbols(med1,med2,circles=0.1,add = TRUE, inches = FALSE,fg = colenv)                                  
        for (i in 1:nrow(matrixdata))
          {
            x <- solve(matrix(c(-med2,med1,med1,med2),nrow=2)
                ,matrix(c(0,med2*coordgenotype[i,dimension2]+med1*coordgenotype[i,dimension1]),ncol=1))
            segments(coordgenotype[i,dimension1],coordgenotype[i,dimension2],x[1],x[2],lty="dotted")
          }                        
      }      
      if (TypeGraph == 8) # Ranking environments
      {
        points(coordgenotype[,dimension1],coordgenotype[,dimension2],col=colgenotype,cex=.8)      
        points(centro[1], centro[2], pch = 18, col = "black")              
        text(coordenviroment[,dimension1],coordenviroment[,dimension2],labels=labelenv,col=colenv,cex=1)
        med1 <<- mean(coordenviroment[,dimension1])
        med2 <<- mean(coordenviroment[,dimension2])
        abline(a=0,b=med2/med1,col=colenv,lty="solid",cex=2)                                                         
        abline(a=0,b=-med1/med2,col=colenv,lty="solid",cex=2)                                                                 
        symbols(med1,med2,circles=0.1,add = TRUE, inches = FALSE,fg = colenv)                  
        mod <<- max((coordenviroment[,dimension1]^2+coordenviroment[,dimension2]^2)^0.5)
        xcoord <<- sign(med1)*(mod^2/(1+ med2^2/med1^2))^0.5
        ycoord <<- (med2/med1)*xcoord
        arrows(centro[1], centro[2], xcoord, ycoord , col = colenv,lty = "solid", length = 0.1)        
        radio <<- ((xcoord-med1)^2+(ycoord-med2)^2)^0.5/3
        for (i in 1:8) symbols(xcoord,ycoord,circles=radio*i,add = TRUE,inches = FALSE,fg = "gray")                 
      }      
      if (TypeGraph == 10) # Ranking genotypes
      {
        med1 <<- mean(coordenviroment[,dimension1])
        med2 <<- mean(coordenviroment[,dimension2])        
        abline(a=0,b=med2/med1,col=colgenotype,lty="solid",cex=2) 
        abline(a=0,b=-med1/med2,col=colgenotype,lty="solid",cex=2)         
        text(coordgenotype[,dimension1],coordgenotype[,dimension2],labels=labelgen,col=colgenotype,cex=1)      
        text(coordenviroment[,dimension1],coordenviroment[,dimension2],labels=labelenv,col=colenv,cex=0.8)        
        coordx <<- 0        
        coordy <<- 0
        for (i in 1:nrow(matrixdata))
          {
            x <- solve(matrix(c(-med2,med1,med1,med2),nrow=2)
                ,matrix(c(0,med2*coordgenotype[i,dimension2]+med1*coordgenotype[i,dimension1]),ncol=1))
            if (sign(x[1])==sign(med1)) 
            {
              if(abs(x[1])> abs(coordx))
              {
                coordx <- x[1]
                coordy <- x[2]
              }
            }
          }  
        arrows(centro[1], centro[2], coordx, coordy , col = colgenotype,lty = "solid", length = 0.1) 
        radio <<- ((coordx-med1)^2+(coordy-med2)^2)^0.5/3
        for (i in 1:10) symbols(coordx,coordy,circles=radio*i,add = TRUE,inches = FALSE,fg = "gray")                                                       
      }            
    } 
    # 
    # Función biplot en 3D
    # 
    Biplot3D <- function()
    {
      dimensions <- 1:3
      rgl.clear("all")              
      rgl.bg(sphere = TRUE, color = c("whitesmoke", "gray90"), lit = FALSE)              
      rgl.light()              
      points3d(coordgenotype[,1],coordgenotype[,2],coordgenotype[,3],col=colgenotype)              
      text3d(coordgenotype[,1], coordgenotype[,2], coordgenotype[,3],labelgen,col=colgenotype,cex=.8)
      text3d(coordenviroment[,1],coordenviroment[,2],coordenviroment[,3],labelenv,col=colenv,cex=1)              
      aspect3d("iso")
      lims <- par3d("bbox")
      segments3d(matrix(c(lims[1], lims[3], lims[5], lims[2], 
                 lims[3], lims[5], lims[1], lims[3], lims[5], lims[1], 
                 lims[4], lims[5], lims[1], lims[3], lims[5], lims[1], 
                 lims[3], lims[6]), byrow = TRUE, ncol = 3), col = "gray60")
      text3d(matrix(c((lims[1] + lims[2])/2, lims[3], lims[5], 
                    lims[1], (lims[3] + lims[4])/2, lims[5], lims[1], 
                    lims[3], (lims[5] + lims[6])/2), byrow = TRUE, nrow = 3), 
                    texts = paste("Dimension ", dimensions), col = "gray60", 
                    family = "sans", font = 1, cex = 1)              
      if (tclvalue(showguidelines)== "1") axes3d()
      for (i in 1:(dim(coordenviroment)[1])) 
      {
        linea <- rbind(coordenviroment[i, ], c(0, 0, 0))
        segments3d(linea[, 1], linea[, 2], linea[, 3], col= colenv)
      }
      if (tclvalue(showtitle) == "1") title3d(wintitle, color = "black", family = "sans", font = 2, cex = 1)
      start <- proc.time()[3]
      while (proc.time()[3] - start < 0.75) 
      {
      }
      start <- proc.time()[3]
      while ((i <- 36 * (proc.time()[3] - start)) < 360) rgl.viewpoint(i, 
              15 - (i - 90)/4, zoom = (if (i < 180)(i + 1)^-0.5 else (360 - i + 1)^-0.5))
              rgl.viewpoint(zoom = 1)
    }       
    #
    # Ventana de selección de genotipos    
    #
    SelectGenotype <- function()
    {
        wingenotype <- tktoplevel()
        tkwm.title(wingenotype,"Select a Genotype")
        combogenotype <- tkwidget(wingenotype,"ComboBox",editable=FALSE,values=levels(labelgen),width=20)
        onOK <- function()
        {
          vgenotype <<- as.numeric(tclvalue(tcl(combogenotype,"getvalue")))+1            
          tkdestroy(wingenotype)
        }
        onCancel <- function()
        {
          vgenotype <<- -1
          tkdestroy(wingenotype)
        }        
        OK.but     <-tkbutton(wingenotype,text="   OK   ",command=onOK)
        Cancel.but <-tkbutton(wingenotype,text=" Cancel ",command=onCancel)
        tkgrid(tklabel(wingenotype,text="    "))                        
        tkgrid(tklabel(wingenotype,text="    "))                        
        tkgrid(tklabel(wingenotype,text="Select a Genotype:    "),combogenotype)            
        tkgrid(tklabel(wingenotype,text="    "))                                          
        tkgrid(OK.but,Cancel.but)
        tkgrid(tklabel(wingenotype,text="    "))
        tkfocus(wingenotype)
        tkwait.window(wingenotype)
    }                               
    #
    # Ventana de selección de ambientes    
    #
    SelectEnvironment <- function()
    {
      winenvironment <- tktoplevel()
      tkwm.title(winenvironment,"Select an Environment")
      comboenvironment <- tkwidget(winenvironment,"ComboBox",editable=FALSE,values=labelenv,width=20)
      onOK <- function()
      {
        venvironment <<- as.numeric(tclvalue(tcl(comboenvironment,"getvalue")))+1            
        tkdestroy(winenvironment)
      }
      onCancel <- function()
      {
        venvironment <<- -1
        tkdestroy(winenvironment)
      }        
      OK.but     <-tkbutton(winenvironment,text="   OK   ",command=onOK)
      Cancel.but <-tkbutton(winenvironment,text=" Cancel ",command=onCancel)
      tkgrid(tklabel(winenvironment,text="    "))                        
      tkgrid(tklabel(winenvironment,text="    "))                        
      tkgrid(tklabel(winenvironment,text="Select an Environment:    "),comboenvironment)            
      tkgrid(tklabel(winenvironment,text="    "))                                          
      tkgrid(OK.but,Cancel.but)
      tkgrid(tklabel(winenvironment,text="    "))
      tkfocus(winenvironment)
      tkwait.window(winenvironment)
    }                    
    #
    # Ventana de selección de ambientes    
    #
    SelectTwoGenotype <- function()
    {
      winEnvGen <- tktoplevel()
      tkwm.title(winEnvGen,"Select Genotypes")
      vgenotype1 <<- -1
      vgenotype2 <<- -1      
      combogenotype1 <- tkwidget(winEnvGen,"ComboBox",editable=FALSE,values=levels(labelgen),width=20)
      combogenotype2 <- tkwidget(winEnvGen,"ComboBox",editable=FALSE,values=levels(labelgen),width=20)      
      onOK <- function()
      {        
        vgenotype1 <<- as.numeric(tclvalue(tcl(combogenotype1,"getvalue")))+1            
        vgenotype2 <<- as.numeric(tclvalue(tcl(combogenotype2,"getvalue")))+1                    
        tkdestroy(winEnvGen)
      }
      onCancel <- function()
      {
        vgenotype1 <<- -1
        vgenotype2 <<- -1
        tkdestroy(winEnvGen)
      }        
      OK.but     <-tkbutton(winEnvGen,text="   OK   ",command=onOK)
      Cancel.but <-tkbutton(winEnvGen,text=" Cancel ",command=onCancel)
      tkgrid(tklabel(winEnvGen,text="Select two genotypes to compare:    "))                  
      tkgrid(tklabel(winEnvGen,text="    "))                        
      tkgrid(tklabel(winEnvGen,text="    "))                        
      tkgrid(tklabel(winEnvGen,text="Genotype 1: "),combogenotype1)
      tkgrid(tklabel(winEnvGen,text="Genotype 2: "),combogenotype2)
      tkgrid(tklabel(winEnvGen,text="    "))                                          
      tkgrid(OK.but,Cancel.but)
      tkgrid(tklabel(winEnvGen,text="    "))
      tkfocus(winEnvGen)
      tkwait.window(winEnvGen)
    }             
    # 
    # Guarda la imagen como un JPG
    #
    SaveFileJPG <- function() 
    {
      FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Jpeg files} {.jpg .jpeg}} {{All files} *}"))
      if (nchar(FileName)) 
        {
          nn <- nchar(FileName)
          if (nn < 5 || substr(FileName, nn - 3, nn) != ".jpg") FileName <- paste(FileName, ".jpg", sep = "")
          jpeg(FileName, width = 8, height = 8, units = "in",restoreConsole = FALSE, res = 96, quality = 50)
          plotFunctiond (screen = FALSE)
          dev.off()
        }
    }  
    #
    # Guarda la imagen como metafile
    # 
    SaveFileMetafile <- function() 
     {
      FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Metafiles} {.wmf}} {{All files} *}"))
      if (nchar(FileName)) 
        {
         nn <- nchar(FileName)
         if (nn < 5 || substr(FileName, nn - 3, nn) != ".wmf") FileName <- paste(FileName, ".wmf", sep = "")
         win.metafile(FileName, width = 8, height = 8, restoreConsole = FALSE)
         plotFunctiond (screen = FALSE)
         dev.off()
        }
    }
    #
    # Guarda la imagen como Postscript
    #
    SaveFilePostscript <- function() 
    {
     FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Postscript files} {.ps}} {{All files} *}"))
     if (nchar(FileName)) 
      {
        nn <- nchar(FileName)
        if (nn < 4 || substr(FileName, nn - 2, nn) != ".ps") FileName <- paste(FileName, ".ps", sep = "")
        postscript(file = FileName, width = 8, height = 8, horizontal = FALSE, onefile = FALSE, paper = "default", 
        family = "URWHelvetica")
        plotFunctiond (screen = FALSE)
        dev.off()
      }
    }                                       
    #
    # Guarda la imagen como PDF
    #
    SaveFilePDF <- function() 
    {
     FileName <- tclvalue(tkgetSaveFile(filetypes = "{{PDF files} {.pdf}} {{All files} *}"))
     if (nchar(FileName)) 
       {
         nn <- nchar(FileName)
         if (nn < 5 || substr(FileName, nn - 3, nn) != ".pdf") FileName <- paste(FileName, ".pdf", sep = "")
         pdf(FileName, width = 7, height = 7)
         plotFunctiond (screen = FALSE)
         dev.off()
       }
    }     
    #
    # Guarda la imagen como bmp
    #
    SaveFileBmp <- function() 
    {
        FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Bitmap files} {.bmp}} {{All files} *}"))
        if (nchar(FileName)) 
        {
         nn <- nchar(FileName)
         if (nn < 5 || substr(FileName, nn - 3, nn) != ".bmp") FileName <- paste(FileName, ".bmp", sep = "")
         bmp(FileName, width = 8, height = 8, units = "in", 
         restoreConsole = FALSE, res = 96)
         plotFunctiond (screen = FALSE)
         dev.off()
      }
    }
    #
    # Guarda la imagen como png
    #
    SaveFilePng <- function() 
    {
     FileName <- tclvalue(tkgetSaveFile(filetypes = "{{Png files} {.png}} {{All files} *}"))
     if (nchar(FileName)) 
       {
         nn <- nchar(FileName)
         if (nn < 5 || substr(FileName, nn - 3, nn) != ".png") FileName <- paste(FileName, ".png", sep = "")
         png(FileName, width = 8, height = 8, units = "in", 
         restoreConsole = FALSE, res = 96)
         plotFunctiond (screen = FALSE)
         dev.off()
       }
    }             
    #
    #  Guarda la imagen como TeX
    #
    SaveFileTeX <- function() 
      {
        FileName <- tclvalue(tkgetSaveFile(filetypes = "{{TeX files} {.tex}} {{All files} *}"))
        if (nchar(FileName)) 
        {
          nn <- nchar(FileName)
          if (nn < 5 || substr(FileName, nn - 3, nn) != ".tex") 
             FileName <- paste(FileName, ".tex", sep = "")
             pictex(FileName, width = 8, height = 8, debug = FALSE, 
             bg = "white", fg = "black")
             plotFunctiond (screen = FALSE)
             dev.off()
        }
    }                         
    #
    # Función que imprime
    #
    Print <- function() 
     {
       try(win.print(), silent = TRUE)
       if (geterrmessage() != "Error in win.print() : unable to start device devWindows\n") 
       {
         plotFunctiond (screen = FALSE)
         dev.off()
       }
    }                                     
    #
    # Función que se ejecuta al pulsar el botón start 
    #
    OnStart <- function()
    {
      tkdestroy(winprincipal)
      wselectfile <- tktoplevel()
      tkwm.title(wselectfile,"GGE Biplot Data")
      OnOpenData <- function()
      {
        tkdestroy(wselectfile)
        file <- tclvalue(tkgetOpenFile())
        if (!length(file))
          return()
        channel <- odbcConnectExcel(file)
        hojastabla <- sqlTables(channel)$TABLE_NAME
        whoja <- tktoplevel()
        tkwm.title(whoja,"Select one table")
        tl <- tklistbox(whoja, height = 4, selectmode = "single",background = "white")
        tkgrid(tklabel(whoja, text = "Select one table:"))
        tkgrid(tl)
        for (k in (1:(dim(sqlTables(channel))[1]))) tkinsert(tl, "end", hojastabla[k])
        tkselection.set(tl, 0)
        OnOKHoja <- function() 
        {
          hojaChoice <<- hojastabla[as.numeric(tkcurselection(tl)) + 1]
          matrizd <<- sqlFetch(channel, hojaChoice) 
          odbcClose(channel)
          tkdestroy(whoja)
        }
        OKHoja.but <- tkbutton(whoja, text = "   OK   ",command = OnOKHoja)
        tkgrid(OKHoja.but)
        tkfocus(whoja) 
        tkwait.window(whoja)     
        mdata <<- matrizd          
        winmodel <- tktoplevel()
        tkwm.title(winmodel, "Model Selection")
        comboscaling <- tkwidget(winmodel,"ComboBox",editable=FALSE,values=datascaling,width=30)
        defaultscaling <- tclVar(optionscaling)
        tkconfigure(comboscaling, textvariable = defaultscaling)
        comboscentering <- tkwidget(winmodel,"ComboBox",editable=FALSE,values=datacentering,width=30)
        defaultcentering <- tclVar(optioncentering)
        tkconfigure(comboscentering, textvariable = defaultcentering)
        comboSVP <- tkwidget(winmodel,"ComboBox",editable=FALSE,values=dataSVP,width=30)
        defaultSVP <- tclVar(optionSVP)
        tkconfigure(comboSVP, textvariable = defaultSVP)       
        OnOKModelSelection <- function()
        {         
          #          
          # Crea la ventana y el menú.
          #
          optioncentering <<- datacentering[as.numeric(tclvalue(tcl(comboscentering,"getvalue")))+1] 
          optionscaling <<- datascaling[as.numeric(tclvalue(tcl(comboscaling,"getvalue")))+1]              
          optionSVP <<- dataSVP[as.numeric(tclvalue(tcl(comboSVP,"getvalue")))+1]                 
          Models()
          tkdestroy(winmodel)
          winplot <- tktoplevel()                
          tkwm.title(winplot, "Graph")          
          img <- tkrplot(winplot, fun = plotFunctiond,hscale = 1.5,vscale = 1.5)
          tkpack(img, expand = "TRUE", fill = "both")                                    
          topMenu <- tkmenu(winplot)
          tkconfigure(winplot,menu=topMenu)
          menuFile <- tkmenu(topMenu,tearoff=FALSE)
          menuView <- tkmenu(topMenu,tearoff=FALSE)          
          menuBiplotTools <- tkmenu(topMenu,tearoff=FALSE)
          menuFormat <- tkmenu(topMenu,tearoff=FALSE)
          menuChangeColor <- tkmenu(topMenu,tearoff=FALSE)
          menuRank <- tkmenu(topMenu,tearoff=FALSE)          
          menuModels <- tkmenu(topMenu,tearoff=FALSE)
          menuBiplot <- tkmenu(topMenu,tearoff=FALSE)          
          menuDividedBy <- tkmenu(topMenu,tearoff=FALSE)
          menuCenteredBy <- tkmenu(topMenu,tearoff=FALSE)          
          menuSVP <- tkmenu(topMenu,tearoff=FALSE) 
          menuSaveAs <- tkmenu(topMenu,tearoff=FALSE)
          tkadd(menuFile,"command",label="Open log file",command=function()  
            {
              Addfile()
            })          
          tkadd(menuFile,"separator")                        
          tkadd(menuFile,"command",label="Copy image",command=function()
            {
              tkrreplot(img)
            })
          tkadd(menuFile,"cascade",label="Save image",menu=menuSaveAs)            
          tkadd(menuSaveAs,"command",label="PDF file",command=function()
            {
              SaveFilePDF()
            })          
          tkadd(menuSaveAs,"command",label="Postscript file",command=function()
            {
              SaveFilePostscript()
            })                                
          tkadd(menuSaveAs,"command",label="Metafile",command=function()
            {
              SaveFileMetafile()
            })          
          tkadd(menuSaveAs,"command",label="Bmp file",command=function()
            {
              SaveFileBmp()
            })           
          tkadd(menuSaveAs,"command",label="Png file",command=function()
            {                       
              SaveFilePng()
            })
          tkadd(menuSaveAs,"command",label="Jpg/Jpeg file",command=function()
            {                       
              SaveFileJPG()
            })            
          tkadd(menuSaveAs,"command",label="TeX file",command=function()
            {                       
              SaveFileTeX()
            })                        
          tkadd(menuFile,"command",label="Print image",command=function()  
            {
              Print()
            })                                
          tkadd(menuFile,"separator")                                                                      
          tkadd(menuFile,"command",label="Exit",command=function()
            {
              tkdestroy(winplot)
            })
          tkadd(menuBiplot,"radiobutton",label="PC1 vs. PC2 (Primary)",variable=vaxis,value = "0",command=function()
            {              
              dimension1 <<- 1
              dimension2 <<- 2
              tkrreplot(img)              
            })            
          tkadd(menuBiplot,"radiobutton",label="PC3 vs. PC4",variable=vaxis,value = "1",command=function()
            {              
              dimension1 <<- 3
              dimension2 <<- 4
              tkrreplot(img)              
            })                        
          tkadd(menuBiplot,"radiobutton",label="PC5 vs. PC6",variable=vaxis,value = "2",command=function()
            {              
              dimension1 <<- 5
              dimension2 <<- 6
              tkrreplot(img)              
            })                                    
          tkadd(menuBiplot,"separator")                                          
          tkadd(menuBiplot,"radiobutton",label="PC1 vs. PC3",variable=vaxis,value = "3",command=function()
            {              
              dimension1 <<- 1
              dimension2 <<- 3
              tkrreplot(img)              
            })                        
          tkadd(menuBiplot,"radiobutton",label="PC2 vs. PC3",variable=vaxis,value = "4",command=function()
            {              
              dimension1 <<- 2
              dimension2 <<- 3
              tkrreplot(img)              
            })                                    
          tkadd(menuBiplot,"separator")                                                                                                        
          tkadd(menuBiplot,"command",label="Biplot 3D",command=function()
            {
              Biplot3D()
            })
          tkadd(menuView,"radiobutton",label="Show Both",variable=showboth,value = "0",command=function()
            {
              tkrreplot(img)
            })                
          tkadd(menuView,"radiobutton",label="Show Genotypes",variable=showboth,value = "1",command=function()          
            {
              tkrreplot(img)
            })                      
          tkadd(menuView,"radiobutton",label="Show Environments",variable=showboth,value = "2",command=function()
            {
              tkrreplot(img)
            })  
          tkadd(menuView,"separator")            
          tkadd(menuDividedBy,"radiobutton",label="0.No scaling",variable=scaling,value = "0",command=function()
            {
              optionscaling <<- "0.No scaling" 
              Models()
              tkrreplot(img)             
            })                                                
          tkadd(menuDividedBy,"radiobutton",label="1.Std Deviation (SD)",variable=scaling,value = "1",command=function()
            {
              optionscaling <<- "1.Std Deviation (SD)" 
              Models()
              tkrreplot(img)                           
            })                                                          
          tkadd(menuCenteredBy,"radiobutton",label="0.No centering",variable=centering,value = "0",command=function()
            {
              optioncentering <<-  "0.No centering"
              Models()
              tkrreplot(img)             
            })                                                            
          tkadd(menuCenteredBy,"radiobutton",label="1.Global-Centered E+G+GE",variable=centering,value = "1",command=function()
            {
              optioncentering <<-  "1.Global-Centered E+G+GE"
              Models()
              tkrreplot(img)             
            })                                                                        
          tkadd(menuCenteredBy,"radiobutton",label="2.Tester-Centered G+GE",variable=centering,value = "2",command=function()
            {
              optioncentering <<-  "2.Tester-Centered G+GE"
              Models()
              tkrreplot(img)             
            })           
          tkadd(menuCenteredBy,"radiobutton",label="3.Double-Centered GE",variable=centering,value = "3",command=function()
            {
              optioncentering <<-  "3.Double-Centered GE"
              Models()
              tkrreplot(img)             
            })                       
          tkadd(menuSVP,"radiobutton",label="JK -(Row Metric Preserving)",variable=svp,value = "0",command=function()
            {
              optionSVP <<- "JK -(Row Metric Preserving)" 
              Models()
              tkrreplot(img)             
            })                                                            
          tkadd(menuSVP,"radiobutton",label="GH -(Column Metric Preserving)",variable=svp,value = "1",command=function()
            {
              optionSVP <<- "GH -(Column Metric Preserving)" 
              Models()
              tkrreplot(img)             
            })                                                                        
          tkadd(menuSVP,"radiobutton",label="HJ -(Dual Metric Preserving)",variable=svp,value = "2",command=function()
            {
              optionSVP <<- "HJ -(Dual Metric Preserving)"
              Models()
              tkrreplot(img)             
            })                                                                                    
          tkadd(menuSVP,"radiobutton",label="SQ - Symmetrical",variable=svp,value = "3",command=function()
            {
              optionSVP <<- "SQ - Symmetrical"
              Models()
              tkrreplot(img)             
            })                                                                                                
          tkadd(menuView,"checkbutton",label="Show/Hide Title",variable=showtitle,command=function()
            {
              if (tclvalue(showtitle)== "1") wintitle <<- "GGE Biplot"    
              if (tclvalue(showtitle)== "0") wintitle <<- NULL
              tkrreplot(img)
            })
          tkadd(menuView,"checkbutton",label="Show/Hide Gidelines",variable=showguidelines,command=function() {tkrreplot(img)})            
          tkadd(menuBiplotTools,"command",label="Examine a Genotype",command=function()
            {
              SelectGenotype()
              if (vgenotype == -1) 
                {
                }
              else
                {
                  wintitle <<- "Examine a Genotype"                                  
                  TypeGraph <<- 3
                  tkentryconfigure(menuView, 2, state = "disabled")                  
                  tkentryconfigure(menuView, 1, state = "disabled")                                    
                  tkrreplot(img)
                }              
            })
          tkadd(menuBiplotTools,"command",label="Examine an Environment",command=function()
            {
              SelectEnvironment()
              if (venvironment == -1) 
                {
                }
              else
                {
                  TypeGraph <<- 2
                  wintitle <<- "Examine an Environment"                  
                  tkentryconfigure(menuView, 2, state = "disabled")                  
                  tkentryconfigure(menuView, 1, state = "disabled")                                    
                  tkrreplot(img)
                }
            })          
          tkadd(menuBiplotTools,"separator")
          tkadd(menuBiplotTools,"command",label="Relation among Environments",command=function()
            {
              TypeGraph <<- 4
              showcircles <<- tclVar("1")
              wintitle <<- "Relationship among environments"
              tkentryconfigure(menuView, 2, state = "disabled")                  
              tkentryconfigure(menuView, 1, state = "disabled")                                    
              tkrreplot(img)                              
            })                                
          tkadd(menuBiplotTools,"separator")                                                  
          tkadd(menuBiplotTools,"command",label="Compare two Genotypes",command=function()
            {
              SelectTwoGenotype()
              TypeGraph <<- 5
              wintitle <<- "Compare two Genotypes"
              tkentryconfigure(menuView, 2, state = "disabled")                  
              tkentryconfigure(menuView, 1, state = "disabled")                                
              tkrreplot(img)
            })                               
          tkadd(menuBiplotTools,"separator")            
          tkadd(menuBiplotTools,"command",label="Which Won Where/What",command=function()
            {
              TypeGraph <<- 6
              tkentryconfigure(menuView, 2, state = "disabled")                  
              tkentryconfigure(menuView, 1, state = "disabled")                                
              tkrreplot(img)            
            })                                            
          tkadd(menuBiplotTools,"command",label="Discrimitiveness vs. representativeness",command=function()
            {      
              TypeGraph <<- 7              
              wintitle <<- "Discrimitiveness vs. representativenss"              
              tkentryconfigure(menuView, 2, state = "disabled")                  
              tkentryconfigure(menuView, 1, state = "disabled")                                
              tkrreplot(img)                           
            })                                                        
          tkadd(menuBiplotTools,"command",label="Mean vs. Stability",command=function()
            {      
              wintitle <<- "Mean vs. Stability"              
              TypeGraph <<- 9              
              tkentryconfigure(menuView, 2, state = "disabled")                  
              tkentryconfigure(menuView, 1, state = "disabled")                                
              tkrreplot(img)                                         
            })                                                                    
          tkadd(menuBiplotTools,"separator")           
          tkadd(menuBiplotTools,"cascade",label="Rank Environment/Genotypes",menu=menuRank)                                                                   
          tkadd(menuRank,"radiobutton",label="with ref.to the 'Ideal' Environment",variable=vrank,value = "1",command=function()
            {
              TypeGraph <<- 8
              wintitle <<- "Ranking Environments"
              tkentryconfigure(menuView, 2, state = "disabled")                  
              tkentryconfigure(menuView, 1, state = "disabled")                                
              tkrreplot(img)                                                     
            })                                                                                            
          tkadd(menuRank,"radiobutton",label="with ref.to the 'Ideal' Genotype",variable=vrank,value = "2",command=function()
            {
              TypeGraph <<- 10
              wintitle <<- "Ranking Genotypes"
              tkentryconfigure(menuView, 2, state = "disabled")                  
              tkentryconfigure(menuView, 1, state = "disabled")                                
              tkrreplot(img)                                                                 
            })                                                                                                        
          tkadd(menuBiplotTools,"separator")                     
          tkadd(menuBiplotTools,"command",label="Back to original data",command=function()
            {
              TypeGraph <<- 1  
              Models()
              showboth <- tclVar("0")
              tkentryconfigure(menuBiplotTools, 0, state = "normal")                                           
              tkentryconfigure(menuView, 2, state = "normal")                  
              tkentryconfigure(menuView, 1, state = "normal")                                
              tkrreplot(img)
            })
          tkadd(menuFormat,"command",label="Plot Title",command=function()
            {
              ReturnVal <<- modalDialog("GGE Biplot","Give your biplot a title:  ","")
              if (ReturnVal=="ID_CANCEL")
              return()
              wintitle <<- ReturnVal
              tkrreplot(img)
              tkfocus(winplot)                
            })          
          tkadd(menuFormat,"separator")
          tkadd(menuChangeColor,"command",label="Background",command=function()
            {
              background <<- ChangeColorv(background)
              tkrreplot(img)            
            })
          tkadd(menuChangeColor,"separator")            
          tkadd(menuChangeColor,"command",label="Genotype labels",command=function()
            {
              colgenotype <<- ChangeColorv(colgenotype)
              tkrreplot(img)            
            })          
          tkadd(menuChangeColor,"command",label="Environment labels",command=function()
            {
              colenv <<- ChangeColorv(colenv)
              tkrreplot(img)            
            }) 
          tkadd(menuChangeColor,"separator")            
          tkadd(menuChangeColor,"command",label="Biplot Title",command=function()
            {
              coltitle <<- ChangeColorv(coltitle)
              tkrreplot(img)                        
            })                                 
          tkadd(menuFormat,"cascade",label="Change Color",menu=menuChangeColor) 
          tkadd(menuModels,"cascade",label="Scaled (divided) by",menu=menuDividedBy)                                                                    
          tkadd(menuModels,"cascade",label="Centered by",menu=menuCenteredBy)
          tkadd(menuModels,"cascade",label="S.V.P.",menu=menuSVP)
          tkadd(topMenu,"cascade",label="File",menu=menuFile)         
          tkadd(topMenu,"cascade",label="View",menu=menuView)                   
          tkadd(topMenu,"cascade",label="Biplot Tools",menu=menuBiplotTools)
          tkadd(topMenu,"cascade",label="Format",menu=menuFormat) 
          tkadd(topMenu,"cascade",label="Models",menu=menuModels)                                                 
          tkadd(topMenu,"cascade",label="Biplot",menu=menuBiplot)                                                           
          tkfocus(winplot)                                    
          if (TypeGraph != "1")
          {
          for (temp1 in 5) tkentryconfigure(menuView, temp1, 
                state = "disabled")
          }
        }
        OK.modelselection <- tkbutton(winmodel, text = "    OK    ",command = OnOKModelSelection)
        tkgrid(tklabel(winmodel, text = "SVP:                                                       "),sticky="w")
        tkgrid(comboSVP)
        tkgrid(tklabel(winmodel, text = "                                                           "),sticky="w")
        tkgrid(tklabel(winmodel, text = "Centered By:                                               "),sticky="w")
        tkgrid(comboscentering)
        tkgrid(tklabel(winmodel, text = "                                                           "),sticky="w")
        tkgrid(tklabel(winmodel, text = "Scaled (Divided) By:                                       "),sticky="w")
        tkgrid(comboscaling)
        tkgrid(tklabel(winmodel, text = "                                                           "),sticky="w")
        tkgrid(OK.modelselection)
        tkfocus(winmodel)
      }   
      SelectData.button <- tkbutton(wselectfile, text = "   Open Data   ", command = OnOpenData)      
      tkgrid(tklabel(wselectfile, text = "    "))            
      tkgrid(tklabel(wselectfile, text = "           Select a table:          ",font=fontFixedWidth))
      tkgrid(tklabel(wselectfile, text = "    "))      
      tkgrid(SelectData.button)
      tkfocus(wselectfile)
    }      
    winprincipal <- tktoplevel()
    tkwm.title(winprincipal, "GGE Biplot")
    Start.button <- tkbutton(winprincipal, text = "   Start   ", command = OnStart)
    fontHeading <- tkfont.create(family = "times", size = 24, weight = "bold", slant = "italic")
    fontFixedWidth <- tkfont.create(family = "courier", size = 12)
    tkgrid(tklabel(winprincipal, text = "               GGE BIPLOT               ",font = fontHeading))
    tkgrid(tklabel(winprincipal, text = "    "))
    tkgrid(Start.button)
    tkfocus(winprincipal)
}

