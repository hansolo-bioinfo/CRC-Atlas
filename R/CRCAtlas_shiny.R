library(shiny)
library(shinyjs)
library(ggplot2)
#library(GSVA)

options(shiny.maxRequestSize = 50*1024^2)

# load data
load("data/lr_search_data_light.RData")
lr_result <- read.csv("data/LR_passed_results_valid.csv", check.names = FALSE)
lr_names <- paste0(lr_result$`ligand(L)`, " -> ", lr_result$`receptor(R)`)
gem_names <- names(top50)

# helper function
plotLR <- function(dose, response, n_para) {
    l <- unlist(strsplit(dose, " -> "))[1]
    r <- unlist(strsplit(dose, " -> "))[2]
    
    dose_name <- dose
    resp_name <- response
    
    if (is.na(r)) {
        x <- npn_norm[[l]]
    } else {
        x <- npn_norm[[l]] * npn_norm[[r]]
    }
    y <- gsva[[response]]
    
    denominator <- max(y)-min(y)
    numerator <- min(y)
    y <- (y-min(y))/(max(y) - min(y))
    d <- data.frame(x, y)
    if (n_para == 2) {
        m <- drc::drm(y ~ x, data = d, fct = drc::LL.2())
        pred <- m$predres[, 1]/10 * denominator + numerator
    } else {
        m <- drc::drm(y ~ x, data = d, fct = drc::LL.4())
        pred <- m$predres[, 1] * denominator + numerator
    }
    sigmoid_rse <- round(sqrt(summary(m)$resVar), digits = 4)
    model_glm <- glm(y ~ x, data = d, family = binomial())
    model_nul <- glm(y ~ 1, data = d, family = binomial())
    stat <- anova(model_nul, model_glm, test = "Chisq")
    p <- round(stat$`Pr(>Chi)`[2], digits = 4)
    
    dd <- data.frame(X = x, Y = gsva[[response]], Pred = pred)
    R <- round(cor(dd$X, dd$Y), digits = 4)
    
    corr <- paste0("Corr=", R)
    rse <- paste0("RSE=", sigmoid_rse)
    pval <- paste0("Pval=", p)
    annotation_text <- paste0(corr, "\n", rse, "\n", pval)
    
    f <- ggplot(dd, aes(x = X, y = Y)) + 
        geom_point(alpha = 2/3, size = 2) + 
        theme_bw() + 
        labs(x = dose_name, y = resp_name) + 
        geom_line(data = dd, aes(x = X, y = Pred, col = "red"), linewidth = 1.5) + 
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 16, face = "bold"),
              legend.position = "none") + 
        annotate("text", x = max(dd$X), y = min(dd$Y), label = annotation_text, 
                 hjust = 1, vjust = 0, size = 4.5)
    f
}

CRCAtlasApp <- function(...) {
    ui <- navbarPage(
        title = "CRC Atlas",
        theme = bslib::bs_theme(4),
        
        # ---- Home ----
        tabPanel(
            title = "Home",
            value = "tab_home",
            icon = icon("house"),
            tags$p(
                "CRC Atlas is an interactive web server that is built on R Shiny for researchers to investigate landscape of gene-coexpressing modules (GEM) in curated single cell RNAseq (scRNAseq) data of colorectal cancer. In CRC Atlas, we have collected 279 samples and 626858 cells from 8 scRNAseq datasets. Please refer to citation for detailed information.", 
                strong("(biorxiv link if applicable)")
            ),
            tags$figure(
                class = "centerFigure",
                tags$img(
                    src = "imgs/Welcome Page.png",
                    width = "1350px",
                    height = "450px",
                    style = "display: block; margin-left: auto; margin-right: auto;"
                ),
                tags$figcaption(strong("Overview of the CRC study."), 
                                em("(a) Merged CRC scRNAseq data from eight cohorts. (b) Illustration of our hierarchical Bayesian model. (c) Output of the our model including the tree structured GEMs, the distribution of GEM over cells, ranked genes associated with each GEM. (d) A framework of discovering the relationship between ligand-receptor (LR) and its downstream GEMs."),
                                style = "font-size: 14px;"
                )
            ),
            tags$hr(
                style = "border-color: lightgrey; 10px;"
            ),
            tags$p(
                "CRC Atlas | © DBMI 2023 | University of Pittsburgh",
                tags$img(
                    src = "imgs/Pitt logo.png",
                    width = "40px",
                    height = "40px"
                ),
                HTML(rep("&nbsp;", 45)),
                "Contact: Han Zhang (haz96@pitt.edu), Dr. Lujia Chen (luc17@pitt.edu)"
            )
        ),
        
        # ---- GEM Information ----
        tabPanel(
            title = "GEM Information",
            value = "tab_stat",
            icon = icon("book"),
            tags$p(
                "Look into the distribution of each gene-coexpressing modules (GEM) and associate top weighted genes by selecting GEM in below dropdown menu. Left panel shows the UMAP, top 20 genes and distribution of the inquired GEM over cell subtypes. Right panel displays the dot plots of all useful GEMs of this major cell type."
            ),
            inputPanel(
                selectInput(
                    inputId = "catype_dropdown",
                    label = "Cell Type",
                    choices = c("T/NK"     = "TNK", 
                                "Myeloid"  = "Myeloid", 
                                "Plasma/B" = "PlasmaB", 
                                "Stromal"  = "Stromal"
                    ),
                    width = "400px"
                ),
                selectInput(
                    inputId = "gem_dropdown",
                    label = "GEM",
                    choices = as.character(c(1:85)),
                    width = "400px"
                )
            ),
            actionButton(
                inputId = "inquire_stat_button",
                "Search"
            ),
            mainPanel(
                fluidRow(
                    column(8, uiOutput("displayImg")),
                    column(4, uiOutput("displaydot"))
                )
            )
        ),
        
        # ---- GEM-LR Correlation ----
        tabPanel(
            title = "GEM-LR Correlation",
            value = "tab_cor",
            icon = icon("chart-line"),
            tags$p(
                "Look into the nonlinear Sigmoid relationship between gene-coexpressing modules (GEM) and ligand -> receptor (LR) pair in our curated CRC Atlas (n=151). On the right panel, you can also find the Top 10 significant GEM-LR pairs for your search, ranked by RSE (residual standard error)."
            ),
            sidebarLayout(
                sidebarPanel(
                    radioButtons(inputId = "gem_lr_choice", 
                                 label   = "Choose one:",
                                 choices = list("GEM", "LR pair"),
                                 inline  = TRUE
                    ),
                    selectInput(
                        inputId = "sigmoid_dropdown",
                        label = "",
                        choices = NULL,
                        width = "300px"
                    ),
                    actionButton(
                        inputId = "inquire_sigmoid_button",
                        "Search"
                    )
                ),
                mainPanel(
                    DT::dataTableOutput("sigmoid_table"),
                    plotOutput("gem_lr_correlation")
                )
            )
        ),
        
        # ---- GEM-LR Conditional Independence ----
        tabPanel(
            title = "GEM-LR Causality",
            value = "tab_ci",
            icon = icon("right-left"),
            tags$p(
                "Conditional independence is an important concept in causal analysis, and in our setting, if a pair of highly correlated GEMs becomes independent when conditioning on a ligand -> receptor (LR) pair, the LR is likely involved in transmitting signal between cells expressing the GEMs. Finding a GEM X-LR-GEM Y triplet leads to the hypothesis depicted in below figure. Herein, we provide the triplet results that pass the conditional independence test.",
                strong("(choose at least one)")
            ),
            mainPanel(
                fluidRow(
                    column(8, tags$figure(
                        tags$img(
                            src = "imgs/conditional_independence.png",
                            width = "800px",
                            height = "280px"
                        ),
                        tags$figcaption(strong("Illustration of conditional independence test."), 
                                        style = "font-size: 14px;"
                        )
                    )),
                    column(4, 
                           selectInput(inputId = "gem_x_drop", 
                                       label   = "Choose GEM-X", 
                                       choices = c("", sort(unique(lr_result$`from_GEM(X)`)))),
                           selectInput(inputId = "gem_lr_drop", 
                                       label   = "Choose Ligand -> Receptor", 
                                       choices = c("", sort(unique(lr_names)))),
                           selectInput(inputId = "gem_y_drop", 
                                       label   = "Choose GEM-Y", 
                                       choices = c("", sort(unique(lr_result$`to_GEM(Y)`))))
                    )
                )
            ),
            mainPanel(
                tags$style(type = "text/css", "#tableContainer { width: 150%; }"),
                div(id = "tableContainer",
                    DT::dataTableOutput("ci_table")
                )
            )
        ),
        
        # ---- Deconvolute ----
        tabPanel(
            title = "Deconvolute",
            value = "tab_deconv",
            icon = icon("magnifying-glass"),
            
            tags$p(
                "To deconvolute the bulk-RNAseq expression data using our curated GEM gene list, you will need to upload the gene expression in <.csv> format of which the rowname is gene symbol and colname is sample ID. The data is better deconvoluted if the input matrix is log2(TPM+1) or log2(FKPM+1)."
            ),
            downloadLink("download_example", "Click here to download an example file"),
            fluidRow(
                column(4, img(src = "imgs/Deconv_Fig.png", 
                              width = "400px", height = "440px")),
                column(8, sidebarPanel(
                    fileInput("file1", "Choose CSV File",
                              accept = c(
                                  "text/csv",
                                  "text/comma-separated-values,text/plain",
                                  ".csv")
                    ),
                    downloadButton("downloadData", "Download Deconvoluted File", 
                                   style = "wide:150%;", disabled = TRUE)
                ))
            )
        )
    )
    
    server <- function(input, output, session) {
        addResourcePath("imgs", "inst/www/")
        
        # ---- GEM Information ----
        observeEvent(input$inquire_stat_button, {
            output$displayImg <- renderUI({
                img(src = paste0("imgs/GEM_Statistics/", input$catype_dropdown, "/",
                                 input$catype_dropdown, "_GEM_", 
                                 input$gem_dropdown, ".png"), 
                    height = "700px", width = "700px")
            })
            
            output$displaydot <- renderUI({
                img(src = paste0("imgs/dotplot/", input$catype_dropdown, ".png"),
                    height = "700px", width = "700px")
            })
        })
        
        # ---- GEM-LR Correlation ----
        observe({
            if (input$gem_lr_choice == "GEM") {
                updateSelectInput(session, "sigmoid_dropdown", 
                                  choices = c("", gem_names))
            } else {
                updateSelectInput(session, "sigmoid_dropdown", 
                                  choices = c("", sort(unique(lr_names))))
            }
        })
        
        observeEvent(input$inquire_sigmoid_button, {
            
            if (input$gem_lr_choice == "GEM") {
                db <- top20_lr[top20_lr$GEM == input$sigmoid_dropdown, ]
            } else {
                db <- top20_gem[top20_gem$`ligand -> receptor` == input$sigmoid_dropdown, ]
            }
            db <- db[order(db$`RSE (Sigmoid)`), ]
            if (nrow(db) >= 10) {
                db <- db[1:10, ]
            }
            
            output$gem_lr_correlation <- renderPlot({
                
                f_lst <- list()
                for (i in 1:10) {
                    f_lst[[i]] <- plotLR(dose     = db$`ligand -> receptor`[i], 
                                         response = db$GEM[i],
                                         n_para   = 2)
                }
                f_lst <- ggpubr::ggarrange(plotlist = f_lst, nrow = 2, ncol = 5, align = "h")
                f_lst
            })
            
            output$sigmoid_table <- DT::renderDataTable({
                DT::datatable(db, options = list(pageLength = 5), rownames = FALSE)
            }, rownames = FALSE)
        })
        
        # ---- GEM-LR Conditional Independence ----
        observe({
            
            updateX <- function(tmp_idx) {
                updateSelectInput(session, "gem_x_drop", 
                                  choices = c("", sort(unique(lr_result$`from_GEM(X)`[tmp_idx]))))
            }
            updateY <- function(tmp_idx) {
                updateSelectInput(session, "gem_y_drop", 
                                  choices = c("", sort(unique(lr_result$`to_GEM(Y)`[tmp_idx]))))
            }
            updateLR <- function(tmp_idx) {
                updateSelectInput(session, "gem_lr_drop", 
                                  choices = c("", sort(unique(lr_names[tmp_idx]))))
            }
            
            empty_x <- input$gem_x_drop == ""
            empty_y <- input$gem_y_drop == ""
            empty_lr <- input$gem_lr_drop == ""
            
            if (empty_x & empty_y & empty_lr) {
                showNotification("Please choose at least one.", type = "message")
                valid_idx <- 1:nrow(lr_result)
                updateX(valid_idx)
                updateY(valid_idx)
                updateLR(valid_idx)
            } else if (empty_x & empty_y & !empty_lr) {
                valid_idx <- which(lr_names == input$gem_lr_drop)
                updateX(valid_idx)
                updateY(valid_idx)
            } else if (empty_x & !empty_y & empty_lr) {
                valid_idx <- which(lr_result$`to_GEM(Y)` == input$gem_y_drop)
                updateX(valid_idx)
                updateLR(valid_idx)
            } else if (!empty_x & empty_y & empty_lr) {
                valid_idx <- which(lr_result$`from_GEM(X)` == input$gem_x_drop)
                updateY(valid_idx)
                updateLR(valid_idx)
            } else if (empty_x & !empty_y & !empty_lr) {
                valid_idx <- which(lr_result$`to_GEM(Y)` == input$gem_y_drop & 
                                       lr_names == input$gem_lr_drop)
                updateX(valid_idx)
            } else if (!empty_x & empty_y & !empty_lr) {
                valid_idx <- which(lr_result$`from_GEM(X)` == input$gem_x_drop & 
                                       lr_names == input$gem_lr_drop)
                updateY(valid_idx)
            } else if (!empty_x & !empty_y & empty_lr) {
                valid_idx <- which(lr_result$`to_GEM(Y)` == input$gem_y_drop & 
                                       lr_result$`from_GEM(X)` == input$gem_x_drop)
                updateLR(valid_idx)
            }
        })
        
        output$ci_table <- DT::renderDataTable({
            if (input$gem_x_drop == "") {
                idx_x <- 1:nrow(lr_result)
            } else {
                idx_x <- which(lr_result$`from_GEM(X)` == input$gem_x_drop)
            }
            if (input$gem_y_drop == "") {
                idx_y <- 1:nrow(lr_result)
            } else {
                idx_y <- which(lr_result$`to_GEM(Y)` == input$gem_y_drop)
            }
            if (input$gem_lr_drop == "") {
                idx_lr <- 1:nrow(lr_result)
            } else {
                idx_lr <- which(lr_names == input$gem_lr_drop)
            }
            idx <- intersect(idx_x, idx_y)
            idx <- intersect(idx, idx_lr)
            
            DT::datatable(lr_result[idx, ], escape=FALSE, 
                          options = list(
                              autoWidth = TRUE,
                              columnDefs = list(list(targets = c(7:14), width = '95px')),
                              scrollX = TRUE,
                              pageLength = 10
                          ))
        }, rownames = FALSE)
        
        # ---- Deconvolute ----
        output$download_example <- downloadHandler(
            filename = function() {
                "expression_example.csv" 
            },
            content = function(file) {
                file.copy("data/expression_example.csv", file)
            }
        )
        
        observeEvent(input$file1, {
            
            if (is.null(input$file1)) {
                return(NULL)
            }
            
            inFile <- input$file1
            data <- read.csv(inFile$datapath, row.names = 1, check.names = FALSE)
            deconv_gem <- GSVA::gsva(as.matrix(data), top50)
            
            enable("downloadData")
            showNotification("File is ready for download", type = "message")
            
            # Enable download button
            output$downloadData <- downloadHandler(
                filename = function() {
                    paste("processed-", inFile$name)
                },
                content = function(file) {
                    write.csv(deconv_gem, file, quote = FALSE)
                }
            )
        })
    }
    
    shinyApp(ui, server, ...)
}
