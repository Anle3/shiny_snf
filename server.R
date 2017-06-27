
 
#load scripts needed for running SNF

options(shiny.maxRequestSize=100*1024^2)
shinyServer(function(input,output){
  fl=reactive({
    validate(
      need(length(input$data)> 1, "Please select at least two data types")
    )
   pprocess_files(features_list[input$data ],input$rm,input$rmv,input$im,input$imv,input$nr)

    })
  #generate similarity matrix
  Wf <- reactive({


    makeWf(K=input$K,tau=input$T,alpha=input$alpha,fl())
    })

  #Vector containing group assignements
  groups <- reactive({cluster_groups(Wf(),input$clusters)})

  #Links for network
  
  output$contents <- renderTable(number_of_clusters(Wf()))

  output$heatmap <- renderPlot({display_heatmap(Wf(),groups())})

  groups2<-reactive({
  groups2=order_group(groups(),descending=T,letter=T)
  groups2=as.data.frame(cbind(names(groups2),groups2))
  groups2$ID=1:nrow(groups2)
  colnames(groups2)=c("name","group","ID")
  return(groups2)
  })
  
  output$networkPlot <-
    renderForceNetwork(

      { 
      if(input$fl_nodes=="topLinks"){
       nw= nw_plot_top(Wf(),groups2(),input$topLinks)
      }else if(input$fl_nodes=="knn"){
          nw=nw_plot_knn(Wf(),groups2(),input$knn)
      }
        
      forceNetwork(
        Nodes = groups2(),
        zoom=TRUE,
        Links = nw,
        Source = "source",
        Target = "target",
        Value = "value",
        NodeID = "name",
        Group = "group",
        opacity=1,
        bounded=TRUE,

        colourScale =JS( 'd3.scale.ordinal().domain("a","b","c","d","e","f","g","h").range(["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd","#8c564b","#e377c2","#7f7f7f","#bcbd22","#17becf"])')
      )
    }
    )



  output$downloadWf <- downloadHandler(
    filename="Similarity_matrix.csv" ,
    content = function(file) {
      write.csv(Wf(), file)
    }
  )

  output$downloadGroups <- downloadHandler(
    filename="Group_assignment.csv",
    conten <- function(file){
      write.table(groups(),file)
    }
  )
  output$downloadNw <- downloadHandler(
    filename="interactions.csv",
    conten <- function(file){
      write.table(nw(),file)
    }
  )

 # output$features<-renderTable({
  #  input$features
   # isolate(select_features(fl(),Wf(),input$clusters))
  #}
  #)

})