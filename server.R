#################################################
#               Factor Analysis                 #
#################################################

library("shiny")
library("nFactors")
library("qgraph")
library("corrplot")

shinyServer(function(input, output) {

  Dataset <- reactive({
    if (is.null(input$file)) { return(NULL) }
    else{
    Dataset <- as.data.frame(read.csv(input$file$datapath ,header=TRUE, sep = ","))
    rownames(Dataset) = Dataset[,1]
    Dataset1 = Dataset[,2:ncol(Dataset)]
    #Dataset = t(Dataset)
    return(Dataset1)
    }
  })

fname <- reactive({
  if(length(strsplit(input$fname,',')[[1]])==0){return(NULL)}
  else{
    return(strsplit(input$fname,',')[[1]])
  }
})

# output$table22 <- renderTable ({ 
#   round(cor(Dataset()),2) 
#                                       })
output$corplot = renderPlot({
  
  my_data = Dataset()
  cor.mat <- round(cor(my_data),2)
  corrplot(cor.mat, 
           type = "upper",    # upper triangular form
           order = "hclust",  # ordered by hclust groups
           tl.col = "black",  # text label color
           tl.srt = 45)  

})    
output$table <- renderTable({ Dataset()   })
  
nS = reactive ({    
  
  if (is.null(input$file)) { return(NULL) }
                    else{
ev = eigen(cor(Dataset(), use = 'pairwise.complete.obs'))  # get eigenvalues
ap = parallel(subject=nrow((Dataset())),var=ncol((Dataset())),rep=100,cent=.05);
nS = nScree(ev$values, aparallel= ap$eigen$qevpea);
}
})

output$fselect <- renderUI({ 
  if (is.null(input$file)) { return(NULL) }
  else{
    
  numericInput("fselect", "Number of Factors:", unlist((nS())[1])[3])
  }
  })

fselect = reactive({
  if (is.null(input$file)) { return(NULL) }
  else{
    
  fselect=input$fselect
  return(fselect)
  }
})

fit = reactive ({ 
  if (is.null(input$file)) { return(NULL) }
  else{
    
  fit = factanal(na.omit(Dataset()), fselect() , scores="Bartlett", rotation="varimax");
  return (fit)
  }
  }) 


output$xaxis <- renderUI({
  
  if (is.null(input$file)) { return(NULL) }
  else {
  if(is.null(fname())){
    n =(fselect())
    list = character(0)
    for (i in 1:n) { 
      temp = paste("Factor",i)
      list = c(list, temp)
    }
    
    selectInput("xaxis", "Choose Factor for plotting on X axis",
                list, selected = "Factor 1")
  }else{
      temp = fname()
       selectInput("xaxis", "Choose Factor for plotting on X axis",
                temp, selected = temp[1])
    
  }
    
    
  }
  
  
})

output$yaxis <- renderUI({
  
  if (is.null(input$file)) { return(NULL) }
  else {
    if(is.null(fname())){
      n =(fselect())
      list = character(0)
      for (i in 1:n) { 
        temp = paste("Factor",i)
        list = c(list, temp)
      }
      list2 = setdiff(list,input$xaxis)
      selectInput("yaxis", "Choose Factor for plotting on Y axis",
                  list2, selected = "Factor 2")
    }else{
      temp = fname()
      temp2 = setdiff(temp,input$xaxis)
      selectInput("yaxis", "Choose Factor for plotting on Y axis",
                  temp2, selected = temp2[1])
      
    }
      }
  
})

f1 = reactive({
  if(is.null(fname())){
    f = input$xaxis
    s <- strsplit(f, "[^[:digit:]]")
    solution <- as.numeric(unlist(s))
    solution <- unique(solution[!is.na(solution)])
    return(solution)
  }else{
index <- match(input$xaxis,fname())
return(index)
}
  
})

f2 = reactive({
  if(is.null(fname())){
    f = input$yaxis
    s <- strsplit(f, "[^[:digit:]]")
    solution <- as.numeric(unlist(s))
    solution <- unique(solution[!is.na(solution)])
    return(solution)
  }else{
    index<-match(input$yaxis,fname()) 
    return(index)
    
   }
  
  
})

output$text1 <- renderText({    
  if (is.null(input$file)) { return(NULL) }
else {
      return(paste("Test of the hypothesis that",(fit())$factors,"factors are sufficient."))}
})

output$text2 <- renderText({
  if (is.null(input$file)) { return(NULL) }
 else{
     return(paste("The chi square statistic is",round((fit())$STATISTIC,3),"on",(fit())$dof," degrees of freedom.")) }
                                   })

output$text3 <- renderText({
  if (is.null(input$file)) { return(NULL) }
 else{
   return(paste("The p-value is",round((fit())$PVAL,3)))
 }
  })
#output$text4 <- renderText({ return(paste("Note - Optimal factors from are:",unlist((nS())[1])[3])) })

output$plot1 = renderPlot({
  if (is.null(input$file)) { return(NULL) }
  else{
    
  plotnScree(nS())
  }
  })


output$plot20 = renderPlot({
  if (is.null(input$file)) { return(NULL) }
  else{
    
a = unclass((fit())$loadings)
grp = NULL
for (i in 1:nrow(a)){
  max = max(abs(a[i,]))
  temp0 =  which(abs(a[i,]) == max)
  temp = temp0[1]
  grp = c(grp,temp)
}
grp = matrix(c(grp,seq(1:length(grp))),,2)
rownames(grp) = colnames(Dataset())

gr = vector("list", length = length(table(grp[,1])))
for (i in 1:length(table(grp[,1]))) {
  l1  = grp[(grp[,1] == as.numeric(names(table(grp[,1])[i]))),2]
  gr[[i]][1:length(l1)] = c(l1)   
}

qgraph(cor(Dataset(), use= 'complete.obs'),layout="spring", groups = gr, labels=names(Dataset()), label.scale=F, label.cex = 1, minimum=input$cutoffcorr)
}
})

output$plot2 = renderPlot({
  if (is.null(input$file)) { return(NULL) }
  else{
    
  a0 = (fit())$loadings
  a1 = a0
  for (i in 1:ncol(a1)){  a1[,i] = a0[,i]*(abs(a0[,i]) > input$cutoff)}
  k2 = f1()
  k3 = f2()
  
  factor.plot = function(a0, a1, k2, k3){
    
    load = a0[((a1[, k2] != 0)|(a1[, k3] != 0)), c(k2, k3)]
    
    par(col="black") #black lines in plots
    
    plot(load,type="p",pch=19,col="red", xlim=c(-1, 1), ylim=c(-1, 1),xlab=input$xaxis,ylab = input$yaxis) # set up plot
    
    abline(h=0);abline(v=0)#draw axes
    
    arrows(0,0, x1=load[,1], y1=load[,2], col="blaCK", lwd=1.5);
    
    text(load,labels = rownames(load),cex=1,pos=1)
    
  } # factor.plot() func ends
  
  factor.plot(a0, a1, k2, k3)
  }
})






output$plot3 = renderPlot({
  if (is.null(input$file)) { return(NULL) }
  else{
    
plot(x=(fit())$scores[,(f1())], y=(fit())$scores[,(f2())], type="p", pch=19, col="red",
    xlab = paste0(input$xaxis), ylab = paste0(input$yaxis))   # added this line in edit

text(x=(fit())$scores[,(f1())],y=(fit())$scores[,(f2())],labels=rownames(Dataset()), pos = 2, col="blue", cex=0.8)

abline(h=0); abline(v=0)
}
})

output$loadings <- renderTable({ 
  if (is.null(input$file)) { return(NULL) } else{
  # rownames((fit())$loadings) = colnames(Dataset())  # edit 2
  b2 <- unclass((fit())$loadings); rownames(b2) <- NULL;  
  b1 <- data.frame(colnames(Dataset()), b2);# [2:ncol(Dataset())];rownames(b1) <- colnames(Dataset())  # edit 2  
  #-------#
  if(is.null(fname())){return(b1)}
  else{
    names(b1)[c(-1)]<-fname()
    return(b1)
  }
  
  #-------#
  return(b1) # unclass((fit())$loadings)
  }
  })

mat = reactive({
  fact = (fit())
# SS.loadings= colSums(fact$loading*fact$loading)
# Proportion.Var = colSums(fact$loading*fact$loading)/dim(fact$loading)[1]
# Cumulative.Var= cumsum(colSums(fact$loading*fact$loading)/dim(fact$loading)[1])
# mat = rbind(SS.loadings,Proportion.Var,Cumulative.Var)
laodings = print(fact$loadings, digits=2, cutoff=.25, sort = TRUE)
# out = list(Stat = mat, loadings=laodings)
# return(laodings)

})

output$mat <- renderPrint({ 
  if (is.null(input$file)) { return(NULL) }
  else{ mat() }
  
  })

# uni = reactive({ 
# a = matrix(fit()$uniqueness,1,)
# colnames(a) = rownames(as.matrix(fit()$uniqueness))
# rownames(a) = "Uniqueness"
# return(a)
#  })

output$uni <- renderTable({ 
  if (is.null(input$file)) { return(NULL) }
  else{ 
    # n = ceiling(length(uni())/3)
    # matrix(uni(), ncol = 1)
    data.frame(Variable = rownames(as.matrix(fit()$uniqueness)), Uniqueness = fit()$uniqueness)
            }
})


output$scores <- renderTable({     
  if (is.null(input$file)) { return(NULL) } else{
      # rownames((fit())$scores) = rownames(Dataset()) # edit 3 i made.
      # b0 <- (fit())$scores;   rownames(b0) <- rownames(Dataset()); 
    
      b0 <- data.frame(rownames(Dataset()), (fit())$scores) # else ends    
      
  if(is.null(fname())){return(b0)}
  else{
    names(b0)[c(-1)]<-fname()
    return(b0)
  }
  }
      # unclass((fit())$scores)
      #                             }
})  
  
# my edits 16-sept-2017 below
  output$downloadDataX <- downloadHandler(
    filename = function() { "Fac_scores.csv" },
    content = function(file) {
      write.csv(data.frame(rownames(Dataset()), (fit())$scores), file, row.names = F)
   				 }
	)

output$downloadData <- downloadHandler(
  filename = function() { "Big_five_Survey_data.csv" },
  content = function(file) {
    write.csv(read.csv("data/Big_five_Survey_data.csv"), file, row.names=F, col.names=F)
  }
)
  
})

