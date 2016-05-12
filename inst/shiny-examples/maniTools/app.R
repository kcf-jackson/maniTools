library(shiny)
library(plotly)
library(maniTools)

dr_demo <- function(sim_data, algor, k, d, kernel = 'rbfdot') {
  if (is.null(algor) | is.null(sim_data) | is.null(k) | is.null(d)) return(NULL)
  if (is.na(algor)) return(NULL)
  algor <- toupper(algor)
  if (algor == "PCA") {
    # PCA (on centered and scaled data)
    run_time <- system.time({
      pca_dr <- sim_data$data %>% center_and_standardise() %>% prcomp()
      proj_data <- sim_data$data %*% pca_dr$rotation[,1:2]
    })
  }

  # MDS
  if (algor == "MDS")
    run_time <- system.time({ proj_data <- cmdscale(dist(sim_data$data), k = d) })

  # Isomap
  if (algor == "ISOMAP")
    run_time <- system.time({ proj_data <- RDRToolbox::Isomap(sim_data$data, dims = d, k = k)$dim2 })

  # LLE
  if (algor == "LLE")
    run_time <- system.time({ proj_data <- LLE2(sim_data$data, dim = d, k = k) })

  if (algor == "DIFFUSIONMAP")
    run_time <- system.time({ proj_data <- diffusionMap::diffuse(dist(sim_data$data), neigen = d)$X })

  # t-SNE
  if (algor == "TSNE")
    run_time <- system.time({ proj_data <- tsne::tsne(sim_data$data, k = d) })

  # KernelPCA
  if (algor == "KPCA")
    run_time <- system.time({ proj_data <- kernlab::kpca(sim_data$data, kernel = kernel, features = d)@pcv })

  # SPE
  if (algor == "SPE")
    run_time <- system.time({ proj_data <- spe::spe(sim_data$data, edim = d)$x })

  # Laplacian Eigenmaps
  if (algor == "LE")
    run_time <- system.time({ proj_data <- Laplacian_Eigenmaps(sim_data$data, k = k, d = d)$eigenvectors })

  # HessianLLE
  if (algor == 'HLLE')
    run_time <- system.time({ proj_data <- Hessian_LLE(sim_data$data, k = k, d = d)$projection })

  # LTSA
  if (algor == 'LTSA')
    run_time <- system.time({ proj_data <- Local_TSA(sim_data$data, k = k, d = d) })

  p1 <- plotly_2D(proj_data, sim_data$colors)
  plot_title <- paste(algor, ". Time taken: ", round(run_time[[1]], 3), "s.", sep = "")
  p1 <- layout(p1, title = plot_title)
  list(p1, run_time)
}

# Define UI for application that draws a histogram
ui <- shinyUI(fluidPage(
   # Application title
   titlePanel("Manifold learning in R"),
   tabsetPanel(
     tabPanel("Demo",
        fluidRow(
          column(6,
             plotlyOutput("plot_3d"),
             plotlyOutput("plot_2d")
          ),
          column(6,
             wellPanel(
               style = "background-color: #ff6666;",
               h3("Manifold"),
               fluidRow(
                 column(4, fileInput("file_input", "Load File", accept = c(
                   'text/csv', 'text/comma-separated-values',
                   'text/tab-separated-values', 'text/plain', '.csv', '.tsv')),
                   checkboxInput('header', 'Header', TRUE)
                 ),
                 column(4, radioButtons('sep', 'Separator',
                           c(Comma=',', Semicolon=';', Tab='\t'),
                           ',')
                 ),
                 column(4, radioButtons('quote', 'Quote',
                           c(None='', 'Double Quote'='"', 'Single Quote'="'"),
                           '"')
                 )
               ),
               textOutput("file_text"),
               hr(),
               fluidRow(
                 column(4,
                        selectInput("data_input", "Examples",
                                    choice = c("Swiss Roll", "Swiss Hole", "Corner Planes",
                                               "Punctured Sphere", "Twin Peaks", "Clusters",
                                               "Toroidal Helix", "Gaussian"), selected = "Swiss Roll")
                 ),
                 column(4, numericInput("num_pts", "#Points", 800)),
                 column(4, uiOutput("ui"))
               )
               #submitButton("Load Example")
             ),

             wellPanel(
               style = "background-color: #00b300;",
               h3("Parameters"),
               fluidRow(
                 column(4, numericInput("d", "Target Dimension d", 2, min = 2, max = 2)),
                 column(4, numericInput("k", "Nearest Neighbors k", 8, min = 1)),
                 column(4, selectInput("kernel", "Kernel",
                               choices = c("rbfdot", #"polydot", "vanilladot", "tanhdot",
                                           "laplacedot", "besseldot", "anovadot"#, "splinedot"
                                           )))
                 #column(3, numericInput("sigma", "Sigma", 10.0)),
                 #column(3, numericInput("alpha", "Alpha", 1.0))
               ),
               textOutput("comment_text")
             ),

             wellPanel(
               style = "background-color: #4d94ff;",
               h3("Algorithm"),
               fluidRow(
                 radioButtons("algor", label = "",
                 choices = list("MDS" = "mds", "PCA" = "pca", "Kernel PCA" = "kpca",
                                "ISOMAP" = "isomap", "Diffusion Map" = "diffusionMap",
                                "Laplacian Eigenmaps" = "le", "Locally Linear Embedding (LLE)" = "lle",
                                "Hessian-LLE (HLLE)" = "hlle", "t-SNE" = "tsne",
                                "Stochastic Proximity Embedding (SPE)" = "spe",
                                "Local Tangent Space Alignment (LTSA)" = "ltsa"),
                 inline = TRUE)
               ),
               textOutput("plot_text")
             )
          )
        )
        ),
     tabPanel("Comparison",
        fluidRow(
          checkboxGroupInput("algor_group", label = h3("Algorithms"),
            choices = list("PCA" = 'pca', "MDS" = 'mds', "ISOMAP" = 'isomap',
                           "LLE" = 'lle', "Diffusion Map" = 'diffusionMap',
                           "t-SNE" = 'tsne', "KPCA" = 'kpca', "SPE" = 'spe',
                           "Laplacian Eigenmaps" = 'le', "HLLE" = 'hlle', "LTSA" = 'ltsa'),
            inline = TRUE),
          actionButton("button", "Update")
        ),
        fluidRow(
          column(4, plotlyOutput("c_plot_1")),
          column(8, uiOutput("plot_first_row"))
        ),
        fluidRow(
          uiOutput("plot_second_row")
        ),
        fluidRow(
          uiOutput("plot_third_row")
        ),
        fluidRow(
          uiOutput("plot_fourth_row")
        )
     )
   )
))





server <- shinyServer(function(input, output) {
  data_from_file <- reactive({
    inFile <- input$file_input
    if (is.null(inFile)) return(NULL)
    sim_data <- read.csv(inFile$datapath, header = input$header,
           sep = input$sep, quote = input$quote)
    x <- sim_data[,1]
    y <- sim_data[,2]
    z <- sim_data[,3]
    if (ncol(sim_data) >= 4) {
      scale <- sim_data[,4]
    } else {
      scale <- z
    }
    return(list(data = cbind(x, y, z), colors = scale))
  })

  DR_data <- reactiveValues(simulation = NULL)
  total_time <- reactiveValues(time_taken = NULL)

  output$ui <- renderUI({
    data_param_label <- switch(input$data_input,
                               "Swiss Roll" = "Height",
                               "Swiss Hole" = "Height",
                               "Corner Planes" = "Angles",
                               "Punctured Sphere" = "Z scale",
                               "Twin Peaks" = "Z scale",
                               "Clusters" = "Number of clusters",
                               "Toroidal Helix" = "Sample rate",
                               "Gaussian" = "Sigma")
    initial_value <- switch(input$data_input,
                            "Swiss Roll" = 1,
                            "Swiss Hole" = 1,
                            "Corner Planes" = 45,
                            "Punctured Sphere" = 1,
                            "Twin Peaks" = 1,
                            "Clusters" = 3,
                            "Toroidal Helix" = 1,
                            "Gaussian" = 1)
    numericInput("data_parameter", data_param_label, value = initial_value)
  })

# First tab ================================================================================
   output$plot_3d <- renderPlotly({
    if (!is.null(data_from_file())) {
      sim_data <- data_from_file()
    } else {
      data_f <- switch(input$data_input,
                       "Swiss Roll" = swiss_roll,
                       "Swiss Hole" = swiss_hole,
                       "Corner Planes" = corner_planes,
                       "Punctured Sphere" = punctured_sphere,
                       "Twin Peaks" = twin_peaks,
                       "Clusters" = clusters_3d,
                       "Toroidal Helix" = toroidal_helix,
                       "Gaussian" = gaussian_random_samples)
      sim_data <- data_f(input$num_pts, input$data_parameter)
      DR_data$simulation <- sim_data
    }
    if (is.null(sim_data$data) | (ncol(sim_data$data) < 3)) {
      plotly_empty()
    } else {
      plotly_3D(sim_data)
    }
  })

   output$plot_2d <- renderPlotly({
     if (!is.null(data_from_file())) {
        sim_data <- data_from_file()
        DR_data$simulation <- sim_data
     }
     if (is.null(DR_data$simulation) | (ncol(DR_data$simulation$data) < 3)) {
        plotly_empty()
     } else {
        res <- dr_demo(DR_data$simulation, algor = input$algor,
                        k = input$k, d = input$d, kernel = input$kernel)
        total_time$time_taken <- res[[2]]
        res[[1]]
     }
   })

   output$file_text <- renderText({"Only plots the first 3 dimensions of the data.
     The 4th dimension is used as colors if available; otherwise, the 3rd dimension is used."})
   output$comment_text <- renderText({"The target dimension is fixed at 2."})
   output$plot_text <- renderPrint({
     cat("Time taken:", total_time$time_taken[[1]], "s. \n")
   })



# Second tab ================================================================================
   output$c_plot_1 <- renderPlotly({
     if (!is.null(data_from_file())) {
       sim_data <- data_from_file()
     } else {
       data_f <- switch(input$data_input,
                        "Swiss Roll" = swiss_roll,
                        "Swiss Hole" = swiss_hole,
                        "Corner Planes" = corner_planes,
                        "Punctured Sphere" = punctured_sphere,
                        "Twin Peaks" = twin_peaks,
                        "Clusters" = clusters_3d,
                        "Toroidal Helix" = toroidal_helix,
                        "Gaussian" = gaussian_random_samples)
       sim_data <- data_f(input$num_pts, input$data_parameter)
       DR_data$simulation <- sim_data
     }
     plotly_3D(sim_data)
   })
   output$plot_first_row <- renderUI({
     plot_output_list <- lapply(1:2, function(i) {
        column(6, plotlyOutput(paste0("c_plot_", i + 1)))
     })
     do.call(tagList, plot_output_list)
   })
   output$plot_second_row <- renderUI({
     plot_output_list <- lapply(3:5, function(i) {
       column(4, plotlyOutput(paste0("c_plot_", i + 1)))
     })
     do.call(tagList, plot_output_list)
   })
   output$plot_third_row <- renderUI({
     plot_output_list <- lapply(6:8, function(i) {
       column(4, plotlyOutput(paste0("c_plot_", i + 1)))
     })
     do.call(tagList, plot_output_list)
   })
   output$plot_fourth_row <- renderUI({
     plot_output_list <- lapply(9:11, function(i) {
       column(4, plotlyOutput(paste0("c_plot_", i + 1)))
     })
     do.call(tagList, plot_output_list)
   })


    observeEvent(input$button, {
       algor_list <- input$algor_group
       for (i in 1:11) {
         local({
           local_i <- i + 1
           output[[paste0("c_plot_", local_i)]] <-
             renderPlotly({
               if ((local_i - 1) %in% seq_along(algor_list)) {
                 if (!is.null(data_from_file())) {
                   DR_data$simulation <- data_from_file()
                 }
                 res <- dr_demo(DR_data$simulation, algor = algor_list[local_i - 1],
                                k = input$k, d = input$d, kernel = input$kernel)
                 total_time$time_taken <- res[[2]]
                 res[[1]]
               } else {
                 plotly_empty()
               }
            })
         })
       }
   })

   # output$c_plot_2 <- renderPlotly({
   #   if (!is.null(data_from_file())) {
   #     tmp_data <- data_from_file()
   #     sim_data <- list(data = as.matrix(tmp_data), colors = tmp_data[,3])
   #     DR_data$simulation <- sim_data
   #   }
   #   algor_list <- input$algor_group
   #   if (length(algor_list) >= 1) {
   #      res <- dr_demo(DR_data$simulation, algor = algor_list[[2]],
   #                     k = input$k, d = input$d, kernel = input$kernel)
   #      total_time$time_taken <- res[[2]]
   #      res[[1]]
   #   }
   # })
#    output$c_plot_3 <- renderPlotly({
#      if (!is.null(data_from_file())) {
#        tmp_data <- data_from_file()
#        sim_data <- list(data = as.matrix(tmp_data), colors = tmp_data[,3])
#        DR_data$simulation <- sim_data
#      }
#      res <- dr_demo(DR_data$simulation, algor = 'mds',
#                     k = input$k, d = input$d, kernel = input$kernel)
#      total_time$time_taken <- res[[2]]
#      res[[1]]
#    })
#    output$c_plot_4 <- renderPlotly({
#      if (!is.null(data_from_file())) {
#        tmp_data <- data_from_file()
#        sim_data <- list(data = as.matrix(tmp_data), colors = tmp_data[,3])
#        DR_data$simulation <- sim_data
#      }
#      res <- dr_demo(DR_data$simulation, algor = 'isomap',
#                     k = input$k, d = input$d, kernel = input$kernel)
#      total_time$time_taken <- res[[2]]
#      res[[1]]
#    })
#    output$c_plot_5 <- renderPlotly({
#      if (!is.null(data_from_file())) {
#        tmp_data <- data_from_file()
#        sim_data <- list(data = as.matrix(tmp_data), colors = tmp_data[,3])
#        DR_data$simulation <- sim_data
#      }
#      res <- dr_demo(DR_data$simulation, algor = 'lle',
#                     k = input$k, d = input$d, kernel = input$kernel)
#      total_time$time_taken <- res[[2]]
#      res[[1]]
#    })
#    output$c_plot_6 <- renderPlotly({
#      if (!is.null(data_from_file())) {
#        tmp_data <- data_from_file()
#        sim_data <- list(data = as.matrix(tmp_data), colors = tmp_data[,3])
#        DR_data$simulation <- sim_data
#      }
#      res <- dr_demo(DR_data$simulation, algor = 'diffusionMap',
#                     k = input$k, d = input$d, kernel = input$kernel)
#      total_time$time_taken <- res[[2]]
#      res[[1]]
#    })
#    output$c_plot_7 <- renderPlotly({
#      if (!is.null(data_from_file())) {
#        tmp_data <- data_from_file()
#        sim_data <- list(data = as.matrix(tmp_data), colors = tmp_data[,3])
#        DR_data$simulation <- sim_data
#      }
#      res <- dr_demo(DR_data$simulation, algor = 'tsne',
#                     k = input$k, d = input$d, kernel = input$kernel)
#      total_time$time_taken <- res[[2]]
#      res[[1]]
#    })
#    output$c_plot_8 <- renderPlotly({
#      if (!is.null(data_from_file())) {
#        tmp_data <- data_from_file()
#        sim_data <- list(data = as.matrix(tmp_data), colors = tmp_data[,3])
#        DR_data$simulation <- sim_data
#      }
#      res <- dr_demo(DR_data$simulation, algor = 'kpca',
#                     k = input$k, d = input$d, kernel = input$kernel)
#      total_time$time_taken <- res[[2]]
#      res[[1]]
#    })
#    output$c_plot_9 <- renderPlotly({
#      if (!is.null(data_from_file())) {
#        tmp_data <- data_from_file()
#        sim_data <- list(data = as.matrix(tmp_data), colors = tmp_data[,3])
#        DR_data$simulation <- sim_data
#      }
#      res <- dr_demo(DR_data$simulation, algor = 'spe',
#                     k = input$k, d = input$d, kernel = input$kernel)
#      total_time$time_taken <- res[[2]]
#      res[[1]]
#    })
#    output$c_plot_10 <- renderPlotly({
#      if (!is.null(data_from_file())) {
#        tmp_data <- data_from_file()
#        sim_data <- list(data = as.matrix(tmp_data), colors = tmp_data[,3])
#        DR_data$simulation <- sim_data
#      }
#      res <- dr_demo(DR_data$simulation, algor = 'le',
#                     k = input$k, d = input$d, kernel = input$kernel)
#      total_time$time_taken <- res[[2]]
#      res[[1]]
#    })
#    output$c_plot_11 <- renderPlotly({
#      if (!is.null(data_from_file())) {
#        tmp_data <- data_from_file()
#        sim_data <- list(data = as.matrix(tmp_data), colors = tmp_data[,3])
#        DR_data$simulation <- sim_data
#      }
#      res <- dr_demo(DR_data$simulation, algor = 'hlle',
#                     k = input$k, d = input$d, kernel = input$kernel)
#      total_time$time_taken <- res[[2]]
#      res[[1]]
#    })
#    output$c_plot_12 <- renderPlotly({
#      if (!is.null(data_from_file())) {
#        tmp_data <- data_from_file()
#        sim_data <- list(data = as.matrix(tmp_data), colors = tmp_data[,3])
#        DR_data$simulation <- sim_data
#      }
#      res <- dr_demo(DR_data$simulation, algor = 'ltsa',
#                     k = input$k, d = input$d, kernel = input$kernel)
#      total_time$time_taken <- res[[2]]
#      res[[1]]
#    })
})

# Run the application
shinyApp(ui = ui, server = server)

