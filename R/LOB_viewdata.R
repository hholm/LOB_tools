LOB_viewdata <- function(LOBpeaklist, rawSpec = NULL){

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
          plotOutput('plot',height = "800px",click = "plot_click"),
          title = "M/Z vs RT",
          br(),
          column(3,
                 h4("Lipid Data Viewer"),
                 checkboxInput('x_auto', 'Auto X-Axis',value = TRUE),
                 sliderInput('rt', 'X-Axis Limits (Retention Time)',
                             min=floor(min(run$peakgroup_rt)),
                             max=ceiling(max(run$peakgroup_rt)),
                             value = c(floor(min(run$peakgroup_rt)),ceiling(max(run$peakgroup_rt))),
                             step=10, round=1),
                 checkboxInput('y_auto', 'Auto Y-Axis',value = TRUE),
                 sliderInput('mz', 'Y-Axis Limits (m/z)',
                             min=floor(min(run$LOBdbase_mz)),
                             max=ceiling(max(run$LOBdbase_mz)),
                             value = c(floor(min(run$LOBdbase_mz)),ceiling(max(run$LOBdbase_mz))),
                             step=10, round=1),

                 checkboxInput('text', 'Display Names'),
                 checkboxInput('oxy', 'Toggle Oxidized Compounds '),
                 checkboxInput('false', 'Toggle False Assignments ')
          ),
          column(3,
                 selectInput('class', 'Select Lipid Class', c("All", as.character(unique(run$species))), multiple = TRUE, selected = "All"),
                 selectInput('color', 'Point Color', c('None','Carbon','Double Bonds', 'Degree Oxidation','Species', 'Lipid Class', 'lpSolve Fitted', 'RF_Window', 'Lipid Class','Final Code')),
                 selectInput('total_carbon', 'Acyl Carbon', c('All',unique(run$FA_total_no_C)), multiple = TRUE, selected = "All"),
                 selectInput('total_db', 'Doubled Bonds', c('All',unique(run$FA_total_no_DB)), multiple = TRUE, selected = "All"),
                 selectInput('deg_ox', 'Degree Oxidation', c('All',unique(run$degree_oxidation)), multiple = TRUE, selected = "All"),
                 selectInput('sizebysample', 'Size by Sample', c("None", colnames(run[13:(length(run)-29)])), multiple = FALSE, selected = "None"
                 ),
                 selectInput('plot_extras', 'Plot Extras', c("RTF_Window", "Labels", "Oxidized Labels"), multiple = TRUE)),
          column(6,
                tableOutput(outputId = "info"),
                actionButton(inputId = "Find_ms2",label = "Check for MS2"),
                numericInput(inputId = "ppm",label = "ppm",value = 5,min = 0,step = 0.1),
                numericInput(inputId = "rtspan",label = "Retention Time Window (s)",value = 200,min = 0,step = 30),
                textOutput(outputId = 'no_sel'),
                tableOutput(outputId = 'ms2_table')
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

      #Set object as NULL to prevent crashing
      run_table <- NULL

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

        # To elimate false assignments if desired
        if(input$false){
          data <- data[data$code!="False_Assignment",]
        }

        # only look at certain acyl lengths
        if("All" %in% input$total_carbon != TRUE){
          data <- data[data$FA_total_no_C == as.numeric(input$total_carbon),]
        }

        # Look at certain DB lengths
        if("All" %in% input$total_db != TRUE){
          data <- data[data$FA_total_no_DB == as.numeric(input$total_db),]
        }

        # Look at certain degrees of oxidation
        if("All" %in% input$deg_ox != TRUE){
          data <- data[data$degree_oxidation == as.numeric(input$deg_ox),]
        }

        # Construct inital plot with limits and points
        g <- ggplot(data = data,mapping = aes(x = peakgroup_rt, y = LOBdbase_mz)) +
          geom_point(aes(size=3)) +
          xlab("Retention Time (sec)") +
          ylab("m/z")

        if(input$x_auto==FALSE){
          g <- g + xlim(c(input$rt[1],input$rt[2]))
        }

        if(input$y_auto==FALSE){
          g <- g + ylim(c(input$mz[1],input$mz[2]))
        }

        #adding labels to each point
        if(!is.null(input$plot_extras)){
        if("Labels" %in% input$plot_extras){
          g <- g +
            geom_text(aes(x = peakgroup_rt, y = LOBdbase_mz,
                          label = (paste0(str_extract(FA_total_no_C, "\\d+"), ":", str_extract(FA_total_no_DB, "\\d+"))),
                          hjust = 1, vjust = 2))
        }}

        #adding oxidation labels to each point
        if(!is.null(input$plot_extras)){
          if("Oxidized Labels" %in% input$plot_extras){
            g <- g +
              geom_text(aes(x = peakgroup_rt, y = LOBdbase_mz,
                            label = (paste0(str_extract(FA_total_no_C, "\\d+"), ":", str_extract(FA_total_no_DB, "\\d+"), ":", str_extract(degree_oxidation, "\\d+"))),
                            hjust = 1, vjust = 2))
          }}

        # change point size by sample
        if (input$sizebysample != "None"){
          g <- g +
            geom_point(data = data, mapping =  aes(x = peakgroup_rt, y = LOBdbase_mz))+
            geom_point(aes(size = !!as.symbol(input$sizebysample)))
        }

        # Add colors for carbon number
        if(input$color=="Carbon"){
          g <- g + geom_point(aes(color=as.character(FA_total_no_C), fill = as.character(FA_total_no_C)), size =3) +
            scale_color_manual(values = palette)

        }

        # Add color for DB number
        if(input$color=="Double Bonds"){
          g <- g + geom_point(aes(color=as.character(FA_total_no_DB), fill = as.character(FA_total_no_DB)),size=3) +
            scale_color_manual(values = palette)
        }

        # Add color for degree oxidation
        if(input$color=="Degree Oxidation"){
          g <- g + geom_point(aes(color=as.character(degree_oxidation), fill = as.character(degree_oxidation)),size=3) +
            scale_color_manual(values = palette)
        }

        # Add color for species
        if(input$color=="Species"){
          g <- g + geom_point(aes(color=species, fill = species),size=3) +
            scale_color_manual(values = palette)
        }

        # Add color for lipid class
        if(input$color=="Lipid Class"){
          g <- g + geom_point(aes(color=lipid_class, fill = lipid_class),size=3) +
            scale_color_manual(values = palette)
        }

        #Add color for final code
        if(input$color=="Final Code"){
          g <- g + geom_point(aes(color=as.character(code), fill = as.character(code)),size=3) +
            scale_color_manual(values = c("LP_Solve_Confirmed"="#66CD00", "10%_rtv"="#66CD00","False_Assignment"="#FF3030", "Red"="#FF3030","RTF_Confirmed"="#2aff00", "ms2v"="#0000FF", "5%_rtv"="#2aff00","LP_Solve_Maybe"="#ff9e44", "Double_Peak?"="#ff9e44", "Double Check"="#ff9e44","Unknown"="#000000", "LP_Solve_Failure"="#B22222", "RTF_Failure"="#B22222"))
        }

        # Add colors for lpSolve solutions
        if(input$color=="lpSolve Fitted"){
          g <- g + geom_point(aes(color=as.character(lpSolve), fill = as.character(lpSolve)),size=3) +
            scale_color_manual(values = c("#E8DA1E","#FF3030","#B6EEA6"))
        }

        if(input$color=="RF_Window"){
          g <- g +
            geom_point(aes(color=as.character(Flag), fill = as.character(Flag)),size=3) +
            scale_color_manual(values = c("LP_Solve_Confirmed"="#66CD00", "10%_rtv"="#66CD00","False_Assignment"="#FF3030", "Red"="#FF3030","RTF_Confirmed"="#2aff00", "ms2v"="#0000FF", "5%_rtv"="#2aff00","LP_Solve_Maybe"="#ff9e44", "Double_Peak?"="#ff9e44", "Double Check"="#ff9e44","Unknown"="#000000", "LP_Solve_Failure"="#B22222", "RTF_Failure"="#B22222"))
        }



        if(!is.null(input$plot_extras)){
        if("RTF_Window" %in% input$plot_extras){
          g <- g +
            geom_point(aes(x = peakgroup_rt/DNPPE_Factor*DBase_DNPPE_RF, y = LOBdbase_mz, color =  as.character(Flag)), shape = 3)+
            geom_errorbarh(aes(xmax = peakgroup_rt/DNPPE_Factor*DBase_DNPPE_RF*1.1, xmin = peakgroup_rt/DNPPE_Factor*DBase_DNPPE_RF*0.9, height = 0.2, y = LOBdbase_mz, color = as.character(Flag)))
        }}


        # Add colors for classes
        if(input$color=="Lipid Class"){
          g <- g + geom_point(aes(color=as.character(species), fill = as.character(species)),size=3) +
            scale_color_manual(values = palette)
        }

        # Add compound names
        if (input$text){
          g <- g + geom_text(aes(label=compound_name),hjust=1,vjust=2,size=3)
        }

        output$info <- renderTable({
          # With ggplot2, no need to tell it what the x and y variables are.
          # threshold: set max distance, in pixels
          # maxpoints: maximum number of rows to return
          # addDist: add column with distance, in pixels
          run_table <- nearPoints(data, input$plot_click, threshold = 20, maxpoints = 1,
                                  addDist = TRUE)
          columns <- c("xcms_peakgroup","compound_name","LOBdbase_mz","peakgroup_rt","Flag","lpSolve","code")
          run[which(run$xcms_peakgroup == run_table$xcms_peakgroup),columns[which(columns %in% colnames(run))]]
        }, digits = 5)

        observeEvent(eventExpr = input$Find_ms2, {
          if (is.null(rawSpec)) {
            output$no_sel <- renderText(paste("No rawSpec objected loaded to read MS2 from."))
          }else{
          run_table <- nearPoints(data, input$plot_click, threshold = 20, maxpoints = 1,
                                  addDist = TRUE)
          withProgress(expr = {
            incProgress(amount = 0.5,message = "Working...")
            if(nrow(run_table)==0){
              output$no_sel <- renderText(
                "Please click a feature on the plot to select it first.")

            }else{
              output$no_sel <- renderText(
                paste("Searching for ms2 data for mass",as.character(run_table$LOBdbase_mz),"... please wait."))
              ms2 <- LOB_findMS2(rawSpec = rawSpec,
                                 mz = run_table$LOBdbase_mz,
                                 rt = run_table$peakgroup_rt,
                                 ppm = input$ppm,
                                 rtspan = input$rtspan)
              output$ms2_table <- renderTable({
                ms2[[1]]
                })
              if (class(ms2[[1]]) == "data.frame") {
              output$no_sel <- renderText(
                paste("Searching for ms2 data for mass",as.character(run_table$LOBdbase_mz),"... Done. First file with most scans:",as.character(names(which(table(ms2[[1]][,'file']) == max(table(ms2[[1]][,'file']))))[1]))
              )}

             # output$no_sel <- renderText("Searching for ms2 data... Done!")
            }})
        }})

        g

      })



      # output$bar <- renderPlot({
      #   tidy_run <- run %>%
      #     gather(key = Sample_ID,
      #            value = Peak_Size,
      #            -species, -compound_name, -LOBdbase_mz,
      #            -peakgroup_rt, -degree_oxidation,
      #            -FA_total_no_C, -FA_total_no_DB, -Flag)
      #   tidy_data <- tidy_run %>%
      #     filter(is.na(FA_total_no_C) != TRUE)
      #   tidy_subset <- tidy_data[tidy_data$species==input$species,]
      #
      #   # further subset for specific flags
      #   if(input$flag != "All"){
      #     tidy_subset <- tidy_subset[tidy_subset$Flag == as.character(input$flag),]
      #   }
      #
      #   # if("All"%in%input$species!=TRUE){
      #   #   tidy_data <- tidy_data[tidy_data$species==input$species,]
      #   # }
      #
      #   # To elimate oxy compounds if desired
      #   # if(input$oxy){
      #   #   tidy_data <- tidy_data[tidy_data$degree_oxidation=="0",]
      #   # }
      #
      #   # Construct inital plot with limits and points
      #   b <- ggplot(data = tidy_subset, aes(x = FA_total_no_C,y = Peak_Size))+
      #     geom_bar(stat = "identity", position = "stack") +
      #     scale_color_manual(values = palette)+# Add colors for carbon number
      #     xlab("as.character(input$xaxis)") +
      #     ylab("Peak Size")+
      #     theme(
      #       axis.text.x = element_blank(),
      #       axis.text.y = element_blank(),
      #       axis.ticks = element_blank())
      #
      #   if(input$xaxis=="Carbon"){
      #     b <- ggplot(data = tidy_subset, aes(FA_total_no_C,y = Peak_Size)) +
      #       geom_bar(stat = "identity", position = "stack") +
      #       scale_color_manual(values = palette)+
      #       xlab(paste(input$xaxis)) +
      #       ylab("Peak Size")# Add colors for carbon number
      #
      #   }
      #
      #   # Add color for DB number
      #   if(input$xaxis=="Double Bonds"){
      #     b <- ggplot(data = tidy_subset, aes(FA_total_no_DB, y = Peak_Size)) +
      #       geom_bar(stat = "identity", position = "stack") +
      #       scale_color_manual(values = palette)+
      #       xlab(paste(input$xaxis)) +
      #       ylab("Peak Size")# Add colors for carbon number
      #   }
      #
      #   # Add colors
      #   if(input$fill_colors=="Double Bonds"){
      #     b <- b + geom_point(aes(fill=FA_total_no_DB)) +
      #       scale_color_manual(values = palette)
      #   }
      #
      #
      #
      #   print(b)
      #
      # })
    }
  )
  runApp(app)
}

