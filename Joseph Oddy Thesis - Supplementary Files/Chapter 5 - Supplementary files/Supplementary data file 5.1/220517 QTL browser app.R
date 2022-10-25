library(shiny)
library(ggplot2)
library(qtl)
library(plyr)
library(rsconnect)

# load in data

H18 <- read.cross("csvs",
                     genfile =  "RXC_rqtl_genotypes.csv",
                     phefile =  "220124 H18 data.csv",
                     genotypes=c("A","B"), crosstype = "dh")

H18 <- jittermap(H18, amount = 1e-6)
H18 <- calc.genoprob(H18, step = 2)

H19 <- read.cross("csvs",
                     genfile =  "RXC_rqtl_genotypes.csv",
                     phefile =  "220124 H19 data.csv",
                     genotypes=c("A","B"), crosstype = "dh")

H19 <- jittermap(H19, amount = 1e-6)
H19 <- calc.genoprob(H19, step = 2)

# Define UI for application

ui <- fluidPage(
  fluidRow(
    column(6, align = "center",
      selectInput(
        inputId = "y1",
        label = "H18 variable",
        choices = c("ln(Alanine)" = "ln.Alanine.", "Arginine" = "Arginine",
                    "Asparagine" = "Asparagine", "Aspartic acid" = "Aspartic_acid",
                    "ln(Glutamic acid)" = "ln.Glutamic_acid.",
                    "ln(Glutamine)" = "ln.Glutamine.", "Glycine" = "Glycine",
                    "-1/Isoleucine" = "X.1.Isoleucine", "Leucine" = "Leucine",
                    "ln(Lysine)" = "ln.Lysine.", "Phenylalanine" = "Phenylalanine",
                    "ln(Serine)" = "ln.Serine.", "Threonine" = "Threonine",
                    "ln(Tyrosine)" = "ln.Tyrosine.", "Area" = "Area",
                    "Length" = "Length", "Width" = "Width", "Diameter" = "Diameter",
                    "ln(KHI)" = "ln.KHI.", "Weight" = "Weight", 
                    "-1/(500-HFN)" = "X.1..500.HFN."),
        selected = "Asparagine")),
    column(6, align = "center",
      selectInput(
        inputId = "y2",
        label = "H19 variable",
        choices = c("ln(Alanine)" = "ln.Alanine.", "Arginine" = "Arginine",
                    "Asparagine" = "Asparagine", "Aspartic acid" = "Aspartic_acid",
                    "ln(Glutamic acid)" = "ln.Glutamic_acid.",
                    "ln(Glutamine)" = "ln.Glutamine.", "Glycine" = "Glycine",
                    "-1/Isoleucine" = "X.1.Isoleucine", "Leucine" = "Leucine",
                    "ln(Lysine)" = "ln.Lysine.", "Phenylalanine" = "Phenylalanine",
                    "ln(Serine)" = "ln.Serine.", "Threonine" = "Threonine",
                    "ln(Tyrosine)" = "ln.Tyrosine.", "Area" = "Area",
                    "Length" = "Length", "Width" = "Width", "Diameter" = "Diameter",
                    "ln(KHI)" = "ln.KHI.", "Weight" = "Weight", 
                    "-1/(500-HFN)" = "X.1..500.HFN."),
        selected = "Asparagine"
      )),
    fluidRow(
      column(12, align = "center", h4("H18 - CIM", align = "center"),
             plotOutput(outputId = "cim1")),
    ),
    fluidRow(
      column(12, align = "center", h4("H18 - Single QTL models", align = "center"),
             tableOutput(outputId = "fulltab1"))
    ),
    fluidRow(
      column(5, align = "center", h4("H18 - Additive QTL model", align = "center"),
             tableOutput(outputId = "H18model")),
      column(7, align = "center", h4("H18 - Dropping terms from additive model", align = "center"),
             tableOutput(outputId = "H18modeleach"))
    ),
    fluidRow(
      column(12, align = "center", h4("H19 - CIM", align = "center"),
             plotOutput(outputId = "cim2")),
    ),
    fluidRow(
      column(12, align = "center", h4("H19 - Single QTL models", align = "center"),
             tableOutput(outputId = "fulltab2"))
    ),
    fluidRow(
      column(5, align = "center", h4("H19 - Additive QTL model", align = "center"),
             tableOutput(outputId = "H19model")),
      column(7, align = "center", h4("H19 - Dropping terms from additive model", align = "center"),
             tableOutput(outputId = "H19modeleach"))
    )
  )
)


# Define server ----------------------------------------------------------------

server <- function(input, output, session) {
  
  # reactive elements ----------------------------------------------------------
  cim1 <- reactive({
    cim(H18, pheno.col = input$y1,
        n.marcovar = nrow(summary((scanone(H18, pheno.col = input$y1)),
                                  threshold = 3)), window = 20, method = "hk")
  })
  
  cim2 <- reactive({
    cim(H19, pheno.col = input$y2,
        n.marcovar = nrow(summary((scanone(H19, pheno.col = input$y2)),
                                  threshold = 3)), window = 20, method = "hk")
  })
  
  cimtab1 <- reactive({
    summary((cim(H18, pheno.col = input$y1,
                n.marcovar = nrow(summary((scanone(H18, pheno.col = input$y1)),
                                          threshold = 3)), window = 20, method = "hk")), threshold = 3)
  })
  
  cimtab2 <- reactive({
    summary((cim(H19, pheno.col = input$y2,
                n.marcovar = nrow(summary((scanone(H19, pheno.col = input$y2)),
                                          threshold = 3)), window = 20, method = "hk")), threshold = 3)
  })
  
  
  # output -------------------------------------------------------------
  output$cim1 <- renderPlot({
    plot(cim1(), bandcol="azure3", alternate.chrid = TRUE, main = input$y1, xlab = "Linkage group", ylab = "LOD score")
    abline(h=3, col="red", lty = 2)
  })
  
  output$cim2 <- renderPlot({
    plot(cim2(), bandcol="azure3", alternate.chrid = TRUE, main = input$y2, xlab = "Linkage group", ylab = "LOD score")
    abline(h=3, col="red", lty = 2)
  })
  
  output$fulltab1 <- renderTable({
    closestmarker1 <- find.marker(H18, chr = cimtab1()$chr,
                                  pos = cimtab1()$pos)
    cimtab1 <- cimtab1()
    cimtab1["Closest marker"] <- closestmarker1
    cimtab1$Peak <- rownames(cimtab1)
    colnames(cimtab1)[1]<-paste("Chromosome")
    colnames(cimtab1)[2]<-paste("cM")
    colnames(cimtab1)[3]<-paste("LOD score")
    cimtab1 <- cimtab1[,c(5,4,1,2,3)]
    
    list1 <- list()
    list1.2 <- list()
    for(i in unique(cimtab1()$chr)){
      output1 <- bayesint(cim1(), chr = i, expandtomarkers = TRUE)
      output1$markers <- rownames(output1)
      chrsubset <- subset(cimtab1(), cimtab1()$chr == i)
      qtl <- makeqtl(H18, chr= i, pos = chrsubset$pos, what=("prob"))
      fitqtl <- fitqtl(H18, pheno.col= input$y1, formula = y ~ Q1, qtl= qtl, method = "hk")
      list1[[i]] <- list(output1)
      list1.2[[i]] <- list(fitqtl$result.full)
    }
    df1 <- ldply (list1, data.frame)
    df1$.id <- NULL
    colnames(df1)[1]<-paste("Chromosome")
    colnames(df1)[2]<-paste("cM")
    colnames(df1)[3]<-paste("LOD score")
    colnames(df1)[4]<-paste("Marker")
    df1 <- df1[,c(4,1,2,3)]
    df1.2 <- ldply(list1.2, data.frame)
    
    cimtab1$lowermarker <- df1$Marker[seq(1, length(df1$Marker), 3)]
    cimtab1$uppermarker <- df1$Marker[seq(3, length(df1$Marker), 3)]
    cimtab1$lowermarkerpos <- df1$cM[seq(1, length(df1$Marker), 3)]
    cimtab1$uppermarkerpos <- df1$cM[seq(3, length(df1$Marker), 3)]
    cimtab1$X.var <- df1.2$X.var[seq(1, length(df1.2$X.var), 3)]
    cimtab1$Pvalue.F. <- df1.2$Pvalue.F.[seq(1, length(df1.2$Pvalue.F.), 3)]
    cimtab1
  })
  
  output$fulltab2 <- renderTable({
    closestmarker2 <- find.marker(H19, chr = cimtab2()$chr,
                                  pos = cimtab2()$pos)
    cimtab2 <- cimtab2()
    cimtab2["Closest marker"] <- closestmarker2
    cimtab2$Peak <- rownames(cimtab2)
    colnames(cimtab2)[1]<-paste("Chromosome")
    colnames(cimtab2)[2]<-paste("cM")
    colnames(cimtab2)[3]<-paste("LOD score")
    cimtab2 <- cimtab2[,c(5,4,1,2,3)]
    
    list2 <- list()
    list2.2 <- list()
    for(i in unique(cimtab2()$chr)){
      output2 <- bayesint(cim2(), chr = i, expandtomarkers = TRUE)
      output2$markers <- rownames(output2)
      chrsubset2 <- subset(cimtab2(), cimtab2()$chr == i)
      qtl2 <- makeqtl(H19, chr= i, pos = chrsubset2$pos, what=("prob"))
      fitqtl2 <- fitqtl(H19, pheno.col= input$y2, formula = y ~ Q1, qtl= qtl2, method = "hk")
      list2[[i]] <- list(output2)
      list2.2[[i]] <- list(fitqtl2$result.full)
    }
    df2 <- ldply (list2, data.frame)
    df2$.id <- NULL
    colnames(df2)[1]<-paste("Chromosome")
    colnames(df2)[2]<-paste("cM")
    colnames(df2)[3]<-paste("LOD score")
    colnames(df2)[4]<-paste("Marker") 
    df2 <- df2[,c(4,1,2,3)]
    df2.2 <- ldply(list2.2, data.frame)
    
    cimtab2$lowermarker <- df2$Marker[seq(1, length(df2$Marker), 3)]
    cimtab2$uppermarker <- df2$Marker[seq(3, length(df2$Marker), 3)]
    cimtab2$lowermarkerpos <- df2$cM[seq(1, length(df2$Marker), 3)]
    cimtab2$uppermarkerpos <- df2$cM[seq(3, length(df2$Marker), 3)]
    cimtab2$X.var <- df2.2$X.var[seq(1, length(df2.2$X.var), 3)]
    cimtab2$Pvalue.F. <- df2.2$Pvalue.F.[seq(1, length(df2.2$Pvalue.F.), 3)]
    cimtab2
  }, rownames = FALSE)
  
  output$H18model <- renderTable({
    H18qtl <- makeqtl(H18, chr=cimtab1()$chr, pos = cimtab1()$pos,what=("prob"))
    H18formula <- as.formula(paste0(" y ~ ", paste(H18qtl$altname, collapse=" + ")))
    H18fitqtl <- fitqtl(H18, pheno.col= input$y1, formula = H18formula, qtl= H18qtl, method = "hk", get.ests = T)
    H18fitqtl$result.full
  }, rownames = TRUE)
  
  output$H18modeleach <- renderTable({
    H18qtl <- makeqtl(H18, chr=cimtab1()$chr, pos = cimtab1()$pos,what=("prob"))
    H18formula <- as.formula(paste0(" y ~ ", paste(H18qtl$altname, collapse=" + ")))
    H18fitqtl <- fitqtl(H18, pheno.col= input$y1, formula = H18formula, qtl= H18qtl, method = "hk", get.ests = T)
    H18fitqtl$result.drop
  }, rownames = TRUE)
  
  output$H19model <- renderTable({
    H19qtl <- makeqtl(H19, chr=cimtab2()$chr, pos = cimtab2()$pos,what=("prob"))
    H19formula <- as.formula(paste0(" y ~ ", paste(H19qtl$altname, collapse=" + ")))
    H19fitqtl <- fitqtl(H19, pheno.col= input$y2, formula = H19formula, qtl= H19qtl, method = "hk", get.ests = T)
    H19fitqtl$result.full
  }, rownames = TRUE)
  
  output$H19modeleach <- renderTable({
    H19qtl <- makeqtl(H19, chr=cimtab2()$chr, pos = cimtab2()$pos,what=("prob"))
    H19formula <- as.formula(paste0(" y ~ ", paste(H19qtl$altname, collapse=" + ")))
    H19fitqtl <- fitqtl(H19, pheno.col= input$y2, formula = H19formula, qtl= H19qtl, method = "hk", get.ests = T)
    H19fitqtl$result.drop
  }, rownames = TRUE)
  
}

# Run the application 
shinyApp(ui = ui, server = server)
