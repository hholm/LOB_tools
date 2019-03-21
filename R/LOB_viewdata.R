LOB_viewdata <- function(LOBpeaklist, RT_Factor_Dbase){

  #Make sure we have our librarys loaded
  library(shiny)
  library(tidyverse)
  library(RColorBrewer)

  #Rename our peak list so we can modify it and keep the complete one
  run <- LOBpeaklist


  # Set up our large color pallete
  palette <-c("#E41A1C","#377EB8","#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
              "#999999", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628",
              "#F781BF", "#999999", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
              "#A65628", "#F781BF", "#999999", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",
              "#FFFF33", "#A65628", "#F781BF", "#999999", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
              "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#E41A1C", "#377EB8", "#4DAF4A",
              "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#E41A1C", "#377EB8",
              "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999", "#E41A1C",
              "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999",
              "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
              "#999999", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628",
              "#F781BF", "#999999", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33",
              "#A65628", "#F781BF", "#999999")

  # Define the app
  app=shinyApp(

    ui = fluidPage(

      title = "Lipid Data Viewer",


      tabsetPanel(
        tabPanel(
          plotOutput('plot',height = "700px",click = "plot_click"),
          title = "M/Z vs RT",
          br(),
          column(3,
                 h4("Lipid Data Viewer"),
                 checkboxInput('x_auto', 'Auto X-Axis'),
                 sliderInput('rt', 'X-Axis Limits (Retention Time)',
                             min=floor(min(run$peakgroup_rt)),
                             max=ceiling(max(run$peakgroup_rt)),
                             value = c(floor(min(run$peakgroup_rt)),ceiling(max(run$peakgroup_rt))),
                             step=10, round=1),
                 checkboxInput('y_auto', 'Auto Y-Axis'),
                 sliderInput('mz', 'Y-Axis Limits (m/z)',
                             min=floor(min(run$LOBdbase_mz)),
                             max=ceiling(max(run$LOBdbase_mz)),
                             value = c(floor(min(run$LOBdbase_mz)),ceiling(max(run$LOBdbase_mz))),
                             step=10, round=1),

                 checkboxInput('text', 'Display Names'),
                 checkboxInput('oxy', 'Toggle Oxidized Compounds ')
          ),
          column(4,
                 selectInput('class', 'Select Lipid Class', c("All", as.character(unique(run$species))), multiple = TRUE, selected = "All"),
                 selectInput('color', 'Point Color', c('None','Carbon','Double Bonds','lpSolve Fitted', 'RF_Window', 'Lipid Class')),
                 selectInput('plot_extras', 'Facet Row', c(None='.',"RTF_Window"),
          )),
          column(5,
                 verbatimTextOutput(outputId = "info",placeholder = TRUE)
                 )),

        tabPanel(
          plotOutput('bar'),
          title = "Distributions",
          column(4,
                 selectInput('species', 'Select Species', c("All",as.character(unique(run$species))),multiple = FALSE),
                 selectInput('xaxis', 'X_Variable', c('None','Carbon','Double Bonds','lpSolve Fitted','Lipid Class'), selected = 'Carbon'),
                 selectInput('fill_colors', 'Fill_Colors', c('Carbon', 'Double Bonds', 'Species'))
          ),
          column(4,
                 selectInput('flag', "Flag?", c("All", as.character(unique(run$Flag))), multiple = TRUE, selected = "All"))

        )
      )
    ),


    # Define server logic to draw our plot
    server = function(input, output) {

      # Will update as varibles change
      output$plot <- renderPlot({

        data <- run #so we dont change our intial data

        # To plot all data
        if("All"%in%input$class!=TRUE){
          data <- data[which(data$species %in% input$class),]
        }

        # To elimate oxy compounds if desired
        if(input$oxy){
          data <- data[data$degree_oxidation=="0",]
        }

        # Construct inital plot with limits and points
        g <- ggplot(data = data,mapping = aes(x = peakgroup_rt, y = LOBdbase_mz)) +
          geom_point(size=3) +
          xlab("Retention Time (sec)") +
          ylab("m/z")

        if(input$x_auto==FALSE){
        g <- g + xlim(c(input$rt[1],input$rt[2]))
        }

        if(input$y_auto==FALSE){
          g <- g + ylim(c(input$mz[1],input$mz[2]))
        }

        # Add colors for carbon number
        if(input$color=="Carbon"){
          g <- g + geom_point(aes(color=as.character(FA_total_no_C))) +
            scale_color_manual(values = palette)

        }

        # Add color for DB number
        if(input$color=="Double Bonds"){
          g <- g + geom_point(aes(color=as.character(FA_total_no_DB)),size=3) +
            scale_color_manual(values = palette)
        }



        # Add colors for lpSolve solutions
        if(input$color=="lpSolve Fitted"){
          g <- g + geom_point(aes(color=as.character(lpSolve)),size=3) +
            scale_color_manual(values = c("#E8DA1E","#FF3030","#B6EEA6"))
        }

        if(input$color=="RF_Window"){
          g <- g +
            geom_point(aes(color=as.character(Flag)),size=3) +
            scale_color_manual(values = c("#00FF00", "#545454","#E8DA1E","#FF3030", "#0000FF"))
        }

        # if(input$plot_extras == "RTF_Window"){
        #   DNPPE_rt <- LOBpeaklist[LOBpeaklist$compound_name=="DNPPE","peakgroup_rt"]
        #   g <- g + geom_point(aes(x =DBase_DNPPE_RF,y = LOBdbase_mz))
        #   #+ geom_errorbarh(aes(xmax = as.numeric(DNPPE_rt*DNPPE_Factor*1.1), xmin = as.numeric(DNPPE_rt*DNPPE_Factor*0.9), color=as.character(Flag)))
        # }

        # Add colors for classes
        if(input$color=="Lipid Class"){
          g <- g + geom_point(aes(color=as.character(species)),size=3) +
            scale_color_manual(values = palette)
        }

        # Add compound names
        if (input$text){
          g <- g + geom_text(aes(label=compound_name),hjust=1,vjust=2,size=3)
        }

        output$info <- renderPrint({
          # With ggplot2, no need to tell it what the x and y variables are.
          # threshold: set max distance, in pixels
          # maxpoints: maximum number of rows to return
          # addDist: add column with distance, in pixels
          nearPoints(data, input$plot_click, threshold = 20, maxpoints = 1,
                     addDist = TRUE)
        })

        g

      })



      output$bar <- renderPlot({
        tidy_run <- run %>%
          gather(key = Sample_ID,
                 value = Peak_Size,
                 -species, -compound_name, -LOBdbase_mz,
                 -peakgroup_rt, -degree_oxidation,
                 -FA_total_no_C, -FA_total_no_DB, -Flag)
        tidy_data <- tidy_run %>%
          filter(is.na(FA_total_no_C) != TRUE)
        tidy_subset <- tidy_data[tidy_data$species==input$species,]

        # further subset for specific flags
        if(input$flag != "All"){
          tidy_subset <- tidy_subset[tidy_subset$Flag == as.character(input$flag),]
        }

        # if("All"%in%input$species!=TRUE){
        #   tidy_data <- tidy_data[tidy_data$species==input$species,]
        # }

        # To elimate oxy compounds if desired
        # if(input$oxy){
        #   tidy_data <- tidy_data[tidy_data$degree_oxidation=="0",]
        # }

        # Construct inital plot with limits and points
        b <- ggplot(data = tidy_subset, aes(x = FA_total_no_C,y = Peak_Size))+
          geom_bar(stat = "identity", position = "stack") +
          scale_color_manual(values = palette)+# Add colors for carbon number
          xlab("as.character(input$xaxis)") +
          ylab("Peak Size")+
          theme(
            axis.text.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks = element_blank())

        if(input$xaxis=="Carbon"){
          b <- ggplot(data = tidy_subset, aes(FA_total_no_C,y = Peak_Size)) +
            geom_bar(stat = "identity", position = "stack") +
            scale_color_manual(values = palette)+
            xlab(paste(input$xaxis)) +
            ylab("Peak Size")# Add colors for carbon number

        }

        # Add color for DB number
        if(input$xaxis=="Double Bonds"){
          b <- ggplot(data = tidy_subset, aes(FA_total_no_DB, y = Peak_Size)) +
            geom_bar(stat = "identity", position = "stack") +
            scale_color_manual(values = palette)+
            xlab(paste(input$xaxis)) +
            ylab("Peak Size")# Add colors for carbon number
        }

        # Add colors
        if(input$fill_colors=="Double Bonds"){
          b <- b + geom_point(aes(fill=FA_total_no_DB)) +
            scale_color_manual(values = palette)
        }



        print(b)

      })
    }
  )
  runApp(app)
}

