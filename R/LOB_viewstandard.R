#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

### Standard Explorer App ###
# Mostly for quick examination of changes in QCs

### Load Packages ###

LOB_viewstandard <- function(object){

  library(shiny)

  library(xcms)

  library(ggplot2)

  # #Create a value for each standard
  rownames <- c("mz","ppm","rtlow","rthigh")
  DNPPE <- c(875.550487, 2.5, 14, 17)
  DGTSd9 <- c(721.66507, 2.5, 15, 19)
  OleicAcidd9 <- c(290.30509, 2.5, 8.26, 12.26)
  ArachidonicAcidd11 <- c(314.30200, 2.5, 7.13,11.13)
  d7MG181 <- c(381.37042, 2.5,7.42,11.42)
  LysoPEd7181 <- c(487.3524, 2.5,6.14,10.14)
  waxesterd5 <- c(505.56839, 2.5,17.24,21.24)
  d7LysoPC181 <- c(529.39935, 2.5,6.05,10.05)
  d7DG15181 <- c(605.58444, 2.5,15.57,19.57)
  d7PE15181 <- c(711.56642, 2.5,14.71,18.71)
  GluCerad5 <- c(733.63488, 2.5,15.11, 18.99)
  PCd715181 <- c(753.61337, 2.5, 14.59,18.59)
  d5PG16181 <- c(771.59064, 2.5, 13.12, 17.12)
  d5TG <- c(857.83285, 2.5, 21.26, 25.26)
  PC320 <- c(734.56943,2.5, 15, 18)
  PG320 <- c(740.54361,2.5, 13, 17)
  PE320 <- c(692.52248,2.5, 15, 18)

  #Make a list for dropdown
  dropdown <- data.frame(DNPPE,
                         DGTSd9,
                         OleicAcidd9,
                         ArachidonicAcidd11,
                         d7MG181,
                         LysoPEd7181,
                         waxesterd5,
                         d7LysoPC181,
                         d7DG15181,
                         d7PE15181,
                         GluCerad5,
                         PCd715181,
                         d5PG16181,
                         d5TG,
                         PC320,
                         PG320,
                         PE320,
                         row.names = rownames)

  ### Set Up the UI ###
  app=shinyApp(
    ui <- shinyUI(fluidPage(

      titlePanel("Standard Explorer"),

      sidebarLayout(

        sidebarPanel(
          #Create a dropdown menu of standards

          selectInput(inputId = "list",
                      label = "Existing Standard",
                      choices = colnames(dropdown)),

          # Input: MZ
          numericInput(inputId = "mz",
                       label = "m/z",
                       value = NULL,
                       min = 300,
                       max = 2000,
                       step = 0.00001,
                       width = NULL),

          # Input: PPM
          numericInput(inputId = "ppm",
                       label = "Mass Tolerance (ppm)",
                       value = NULL,
                       min = 1,
                       max = 100,
                       step = 1,
                       width = NULL),

          # Input: RT Low
          numericInput(inputId = "rtmin",
                       label = "Retention Time - Min",
                       value = NULL,
                       min = 0,
                       max = 30,
                       step = 0.01,
                       width = NULL),

          # Input: RT High
          numericInput(inputId = "rtmax",
                       label = "Retention Time - Max",
                       value = NULL,
                       min = 0,
                       max = 30,
                       step = 0.01,
                       width = NULL),

          # Input: Samples to Graph
          selectInput(inputId = "graph_num",
                       choices = phenoData(object)@data,
                       multiple = TRUE,
                       label = "Sample Numbers for Chromatogram",
                      ),

          #checkboxInput(inputId = "XYZ", label = "Group peaks together by sample?"),

          #GOOOOOOO!
          actionButton(inputId = "runtest",
                       "Search for Peak IDs"),
          verbatimTextOutput("parameters"),
          actionButton(inputId = "rungraph",
                       "Generate Chromatograms")

        ),

        #Add click button here.

        mainPanel(
          tabsetPanel(
            tabPanel("Peak Table",
                     tableOutput("table"),
                     textOutput("nopeaks")),
            tabPanel("Statistic Plots",
                     verbatimTextOutput("statistics"),
                     plotOutput("rt_graph"),
                     plotOutput("intensity_graph"),
                     plotOutput("mz_graph")),
            tabPanel("Peak",
                     plotOutput("plot2",height = "800px",click = "plot_click",hover = hoverOpts("plot_hover", delay = 100, delayType = "debounce")),
                     plotOutput("plot",height = "300px"),
                     uiOutput("hover_info")
                     )
          )
        )
      )
    )),

    ### server - The code behind the UI
    server <- function(input, output) {

      ### Make the table and Graph ###

      observeEvent(eventExpr = input$rungraph, {
        if(is.null(mz)){
          output$nopeaks <- renderText(
            "Please search for Peak IDs first.")
        }else{

          #Create the graph
        cent_file <- filterFile(object,file = input$graph_num)
        chroms <- plotChromPeakDensity(object = object, mz = c(mzlow,mzhigh),rt = c(seclow, sechigh))

        output$plot <- renderPlot(plot(chroms))
        }

      }
    )

      observeEvent(eventExpr = input$runtest, {

        output$nopeaks <- renderText("")
        output$table <- renderTable("")



        mz <- if(is.na(input$mz) == TRUE){
          dropdown["mz",as.character(input$list)]
        } else {
          input$mz
        }
        ppm <- if(is.na(input$ppm) == TRUE){
          dropdown["ppm",as.character(input$list)]
        } else {
          input$ppm
        }
        rtlow <- if(is.na(input$rtmin) == TRUE){
          dropdown["rtlow",as.character(input$list)]
        } else {
          input$rtmin
        }
        rthigh <- if(is.na(input$rtmax) == TRUE){
          dropdown["rthigh",as.character(input$list)]
        } else {
          input$rtmax
        }

        #Create the text box for parameters
        output$parameters <- renderText(c('Current Settings','\nm/z =',mz,
                                          '\nppm =',ppm,
                                          '\nrtlow =',rtlow,
                                          '\nrthigh =',rthigh))

        #turn our minutes into seconds
        seclow <- (rtlow*60)
        sechigh <- (rthigh*60)

        #create our m/z range based on the ppm
        mzrange <- mz*(0.000001*ppm)
        mzlow <- (mz-mzrange)
        mzhigh <- (mz+mzrange)

        #make a data frame of our sample names

        samplenames <- sampleNames(object)

        samplenamesframe <- data.frame(samplenames,samplenumber =
                                         seq(from=1, to=length(mzXMLfiles)))

        #create + extract a lists of peaks that fit our parameters
        peaks <- chromPeaks(object = object,
                            mz = c(mzlow, mzhigh),
                            rt = c(seclow, sechigh))
        if (hasFeatures(object)) {
          fts <- featureDefinitions(object,
                                    mz = c(mzlow,mzhigh),
                                    rt = c(seclow,sechigh))
        }

        #turn our matrix into a dataframe
        peaksframe <- as.data.frame(peaks)

        #pull out the columns we want
        peaksnumber <- peaksframe[["sample"]]
        peaksmz <- peaksframe[["mz"]]
        peaksrt <- peaksframe[["rt"]]
        peaksintensity <-peaksframe[["into"]]

        #make them into another dataframe
        samplevalues <- data.frame(name = peaksnumber,
                                   mz = peaksmz,
                                   rt = peaksrt,
                                   intensity = peaksintensity)

        #Add the sample names back in
        merged <- merge(samplevalues, samplenamesframe, by.x="name", by.y= "samplenumber")

        #Reorder our coulmns so sample name comes seconds
        reordered <- merged[c(1,5,2,3,4)]

        if(is.na(reordered[1,"name"])== TRUE){
          output$nopeaks <- renderText(
            "No peaks found in object for current settings.")
        }else{

          #Make everything a character so we can add a page break in <- made switch but dont like it

          #if(input$XYZ == TRUE){
          ascharacters <- as.data.frame(lapply(reordered, as.character), stringsAsFactors = FALSE)

          #Done <- head(do.call(rbind, by(ascharacters, reordered$name, rbind, "")), -1 )
          # }else{

          # Done <- reordered
          # }

          Done <- ascharacters

          colnames(Done)[1] <- "Sample Number"
          colnames(Done)[2] <- "Sample Name"
          colnames(Done)[3] <- "m/z"
          colnames(Done)[4] <- "Retention Time"
          colnames(Done)[5] <- "Intensity"

          #Sci format
          Done[5] <- format(as.numeric(Done[[5]]), scientific = TRUE)
          Done[4] <- (as.numeric(Done[[4]])/60)


          output$table <- renderTable(print(Done),
                                      striped = FALSE,
                                      align = 'l',
                                      width = 400,
                                      digits = 5
                                     )
        }
        Extra_line <- NULL
        Extra_points <- NULL
        Extra <- rep("No Group",nrow(peaks))
        if (hasFeatures(object) & nrow(fts)>0) {
          pnames<- gsub(x = rownames(peaks),pattern = c("CP000000"),replacement = "")
          pnames<- gsub(x = pnames,pattern = c("CP00000"),replacement = "")
          pnames<- gsub(x = pnames,pattern = c("CP0000"),replacement = "")
          pnames<- gsub(x = pnames,pattern = c("CP000"),replacement = "")
          pnames<- gsub(x = pnames,pattern = c("CP00"),replacement = "")
          pnames<- gsub(x = pnames,pattern = c("CP0"),replacement = "")
          pnames<- gsub(x = pnames,pattern = c("CP"),replacement = "")
          for (i in length(fts@listData[["peakidx"]])) {
            feat <- fts@listData[["peakidx"]][[i]]
            which <- pnames %in% feat
            Extra[which(pnames %in% feat)] <- fts@rownames[i]
          }
          Extra_line <- aes(group = Extra)
          Extra_points <- aes(color = Extra)
          }else{
          Extra_line <- NULL
          Extra_points <- NULL
        }

        #Create STD deveation stats
        output$rt_graph <- renderPlot(
          ggplot(data = reordered,aes(x = name ,y = rt,group=1)) +
            geom_point(mapping = Extra_points) +
            geom_line(mapping = Extra_line,linetype="dotted") +
            ylab("Retention Time (Seconds)") +
            xlab("Sample Number")
        )
        output$intensity_graph <- renderPlot(
          ggplot(data = reordered,aes(x = name ,y = as.numeric(intensity),group=1)) +
            geom_point(mapping = Extra_points) +
            geom_line(mapping = Extra_line,linetype="dotted") +
            ylab("Intensity") +
            xlab("Sample Number") +
            scale_y_continuous(trans='log10')

        )

        mz_ppm_diff <- reordered$mz-mz
        mz_ppm_diff <- mz_ppm_diff/(0.000001*mz)

        output$mz_graph <- renderPlot(
          ggplot(data = reordered,aes(x = name ,y = mz_ppm_diff,group=1)) +
            geom_point(mapping = Extra_points) +
            geom_line(mapping = Extra_line,linetype="dotted") +
            ylab("Mass / Change (m/z)") +
            xlab("Sample Number")
        )

        output$statistics <- renderText(c('Statistics','\nIntensity Deviation =',sd(as.numeric(reordered$intensity)),
                                          '\nppm =',ppm,
                                          '\nrtlow =',rtlow,
                                          '\nrthigh =',rthigh))

        #Plot grouping info
        cols <- RColorBrewer::brewer.pal(8,name = "Dark2")
        pall <- colorRampPalette(cols)
        pks <- data.frame(peaks)
        pks$feature_group <- Extra
        colors <-pall(nrow(pks))

        output$plot <- renderPlot(
          plotChromPeakDensity(object = object,
                               mz = c(mzlow,mzhigh),
                               rt = c(seclow, sechigh),
                               col = colors,pch = 16)
                                 )
       g <- ggplot() +
          geom_point(data = pks,size = 3,aes(x = rt,y = sample,color = sampleNames(object)[pks$sample]))+
          geom_density(aes(x = pks$rt,y = ..scaled..*length(pks$sample))) +
          theme_minimal() +
          theme(legend.title = element_blank(),legend.position ="none")

        if (hasFeatures(object)) {
          g <- g + geom_rect(fill=alpha("grey",0),alpha = 0.5,aes(xmin=fts$rtmin,xmax=fts$rtmax,ymin=min(pks$sample),ymax=max(pks$sample)))
        }
       cur_x<- ggplot_build(g)$layout$panel_scales_x[[1]]$range$range
       g <- g + xlim(cur_x[1]-20,cur_x[2]+20)
        output$plot2 <- renderPlot(
          g
        )
        output$hover_info <- renderUI({
          hover <- input$plot_hover
          point <- nearPoints(pks, hover, threshold = 5, maxpoints = 1, addDist = TRUE)
          if (nrow(point) == 0) return(NULL)

          # calculate point position INSIDE the image as percent of total dimensions
          # from left (horizontal) and from top (vertical)
          left_pct <- (hover$x - hover$domain$left) / (hover$domain$right - hover$domain$left)
          top_pct <- (hover$domain$top - hover$y) / (hover$domain$top - hover$domain$bottom)

          # calculate distance from left and bottom side of the picture in pixels
          left_px <- hover$range$left + left_pct * (hover$range$right - hover$range$left)
          top_px <- hover$range$top + top_pct * (hover$range$bottom - hover$range$top)

          # create style property fot tooltip
          # background color is set so tooltip is a bit transparent
          # z-index is set so we are sure are tooltip will be on top
          style <- paste0("position:absolute; z-index:100; background-color: rgba(245, 245, 245, 0.85); ",
                          "left:", left_px+2, "px; top:", top_px+50, "px;")

          # actual tooltip created as wellPanel
          wellPanel(
            style = style,
            p(HTML(paste0("<b> Peak: </b>", rownames(point), "<br/>",
                          "<b> MZ: </b>", point$mz, "<br/>",
                          "<b> RT: </b>", point$rt, "<br/>",
                          "<b> FT: </b>", point$feature_group, "<br/>",
                          "<b> Samp: </b>", sampleNames(object)[point$sample], "<br/>")))
          )
        })

        output$info <- renderTable({
          # With ggplot2, no need to tell it what the x and y variables are.
          # threshold: set max distance, in pixels
          # maxpoints: maximum number of rows to return
          # addDist: add column with distance, in pixels
          run_table <-nearPoints(pks, input$plot_click, threshold = 20, maxpoints = 1,
                                 addDist = TRUE)
        }, digits = 5)

      })
    }

    # Run the application
    #shinyApp(ui = ui, server = server)
  )
  runApp(app)
}


