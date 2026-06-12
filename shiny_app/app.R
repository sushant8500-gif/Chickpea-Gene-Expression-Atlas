# ============================================================
# Chickpea Gene Expression Atlas + GWAS Candidate Gene Finder
# ============================================================

setwd("C:/Users/lw267/OneDrive - Tennessee State University/Desktop/Chapter 2_analysis agaiin/Chickpea expression atlas_RK")

library(shiny)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotly)
library(DT)
library(xml2)
library(htmltools)
library(scales)
library(bslib)
library(rsvg)
library(magick)
library(openxlsx)

# ============================================================
# 1. File paths
# ============================================================

expr_file  <- "data/FPKM_File_RK.xlsx"
svg_file   <- "www/Chickpea_gene_expression_atlas_RK.svg"
annot_file <- "data/gene_annotation.xlsx"

if (!file.exists(expr_file)) stop("Expression file not found: data/FPKM_File_RK.xlsx")
if (!file.exists(svg_file)) stop("SVG file not found: www/Chickpea_gene_expression_atlas_RK.svg")
if (!file.exists(annot_file)) stop("Annotation file not found: data/gene_annotation.xlsx")

# ============================================================
# 2. Tissue IDs
# ============================================================

tissues <- c(
  "Androecium", "Bracteole", "Embryo", "Endosperm",
  "Flower_1", "Flower_2", "Flower_3", "Flower_4", "Flower_5",
  "Flower_bud_1", "Flower_bud_2", "Flower_bud_3", "Flower_bud_4",
  "Germinating_Seedling", "Gynoecium", "Mature_leaf", "Nodule",
  "Pedicel", "Petal", "Pod_Shell", "Root", "Root_Hair", "Root_tip",
  "SAM", "Seed_10_dap", "Seed_20_dap", "Seed_30_dap", "Seed_5_dap",
  "Seed_Coat", "Sepal", "Shoot", "Young_leaf"
)

tissue_groups <- list(
  "All tissues" = tissues,
  "Root/Nodule" = c("Root", "Root_Hair", "Root_tip", "Nodule"),
  "Leaf/Shoot/SAM" = c("Mature_leaf", "Young_leaf", "Shoot", "SAM"),
  "Flower organs" = c(
    "Flower_1", "Flower_2", "Flower_3", "Flower_4", "Flower_5",
    "Flower_bud_1", "Flower_bud_2", "Flower_bud_3", "Flower_bud_4",
    "Androecium", "Gynoecium", "Petal", "Sepal", "Pedicel", "Bracteole"
  ),
  "Seed/Pod" = c(
    "Seed_5_dap", "Seed_10_dap", "Seed_20_dap", "Seed_30_dap",
    "Seed_Coat", "Embryo", "Endosperm", "Pod_Shell"
  )
)

# ============================================================
# 3. Read expression matrix
# ============================================================

expr <- read_excel(expr_file)

gene_col <- "Gene_ID"

if (!gene_col %in% colnames(expr)) {
  stop("Gene_ID column not found in FPKM_File_RK.xlsx. Change gene_col in app.R.")
}

missing_excel <- setdiff(tissues, colnames(expr))

if (length(missing_excel) > 0) {
  stop(paste("These tissues are missing in Excel:", paste(missing_excel, collapse = ", ")))
}

expr <- expr %>%
  mutate(across(all_of(tissues), ~ as.numeric(.)))

gene_choices <- sort(unique(expr[[gene_col]]))

# ============================================================
# 4. Read gene annotation file
# ============================================================

gene_annot <- read_excel(annot_file)

required_annot_cols <- c(
  "Gene_ID", "Gene_raw", "Chr", "Start", "End",
  "Strand", "Gene_length_bp", "Annotation", "Dbxref"
)

missing_annot_cols <- setdiff(required_annot_cols, colnames(gene_annot))

if (length(missing_annot_cols) > 0) {
  stop(paste("These columns are missing in gene_annotation.xlsx:", paste(missing_annot_cols, collapse = ", ")))
}

gene_annot <- gene_annot %>%
  mutate(
    Chr = as.character(Chr),
    Start = as.numeric(Start),
    End = as.numeric(End),
    Gene_length_bp = as.numeric(Gene_length_bp)
  )

# ============================================================
# 5. Helper functions
# ============================================================

palette_choices <- list(
  "Yellow-orange-red" = c("#ffffcc", "#ffeda0", "#feb24c", "#f03b20", "#bd0026"),
  "White-yellow-red" = c("#f7fbff", "#ffffb2", "#fd8d3c", "#bd0026"),
  "White-orange-red" = c("#fff7ec", "#fdd49e", "#fc8d59", "#b30000"),
  "Blue-white-red" = c("#2166ac", "#f7f7f7", "#b2182b"),
  "Purple-yellow" = c("#2d004b", "#762a83", "#f7f7f7", "#fdb863", "#e66101"),
  "Viridis" = c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"),
  "Magma" = c("#000004", "#3b0f70", "#8c2981", "#de4968", "#fe9f6d", "#fcfdbf"),
  "Plasma" = c("#0d0887", "#6a00a8", "#b12a90", "#e16462", "#fca636", "#f0f921"),
  "Cividis" = c("#00204c", "#31446b", "#666970", "#958f78", "#c6ba7c", "#ffea46"),
  "Green-yellow-red" = c("#006837", "#78c679", "#ffffbf", "#fdae61", "#a50026"),
  "Light-blue-dark-blue" = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b")
)

expression_to_color <- function(x, max_value, pal_colors) {
  x <- as.numeric(x)
  x[is.na(x)] <- 0
  
  x <- pmin(x, max_value)
  scaled <- x / max_value
  
  pal <- colorRampPalette(pal_colors)(100)
  pal[pmax(1, round(scaled * 99) + 1)]
}

color_svg <- function(gene_row, max_fpkm, pal_colors,
                      stroke_mode = "No stroke",
                      stroke_color = "#333333",
                      stroke_width = 1) {
  
  svg <- read_xml(svg_file)
  
  for (tissue in tissues) {
    
    value <- as.numeric(gene_row[[tissue]])
    color <- expression_to_color(value, max_fpkm, pal_colors)
    
    node <- xml_find_first(svg, paste0("//*[@id='", tissue, "']"))
    
    if (!is.na(xml_name(node))) {
      
      if (stroke_mode == "No stroke") {
        xml_set_attr(
          node,
          "style",
          paste0("fill:", color, ";stroke:none;")
        )
      } else {
        xml_set_attr(
          node,
          "style",
          paste0(
            "fill:", color,
            ";stroke:", stroke_color,
            ";stroke-width:", stroke_width,
            ";"
          )
        )
      }
    }
  }
  
  as.character(svg)
}

parse_gene_input <- function(selected_genes, pasted_genes) {
  pasted <- unlist(strsplit(pasted_genes, "\\s+|,|;"))
  pasted <- pasted[pasted != ""]
  genes <- unique(c(selected_genes, pasted))
  genes <- genes[genes %in% expr[[gene_col]]]
  genes
}

normalize_chr <- function(x) {
  x <- as.character(x)
  x <- gsub("chr", "", x, ignore.case = TRUE)
  x <- gsub("^0+", "", x)
  x
}

# ============================================================
# Fixed 4500 x 3250 atlas export function
# ============================================================

save_atlas_with_legend <- function(svg_text, file, format = "png",
                                   final_width = 4500,
                                   final_height = 3250,
                                   max_fpkm = 100,
                                   pal_colors = c("#f7fbff", "#ffffb2", "#fd8d3c", "#bd0026")) {
  
  tmp_svg <- tempfile(fileext = ".svg")
  tmp_png <- tempfile(fileext = ".png")
  legend_file <- tempfile(fileext = ".png")
  
  writeLines(svg_text, tmp_svg, useBytes = TRUE)
  
  rsvg::rsvg_png(
    svg = tmp_svg,
    file = tmp_png,
    width = final_width
  )
  
  atlas_img <- magick::image_read(tmp_png)
  atlas_img <- magick::image_background(atlas_img, color = "white", flatten = TRUE)
  
  png(
    filename = legend_file,
    width = 1400,
    height = 260,
    res = 150,
    bg = "white"
  )
  
  par(mar = c(4.5, 4, 1, 1), bg = "white")
  
  plot(
    NA,
    xlim = c(0, max_fpkm),
    ylim = c(0, 1),
    axes = FALSE,
    xlab = "",
    ylab = "",
    xaxs = "i",
    yaxs = "i"
  )
  
  pal <- colorRampPalette(pal_colors)(100)
  breaks <- seq(0, max_fpkm, length.out = 101)
  
  for (i in 1:100) {
    rect(
      xleft = breaks[i],
      ybottom = 0,
      xright = breaks[i + 1],
      ytop = 1,
      col = pal[i],
      border = NA
    )
  }
  
  box()
  axis(1, at = pretty(c(0, max_fpkm), n = 5), cex.axis = 0.9)
  mtext("FPKM expression", side = 1, line = 2.7, cex = 1.05)
  
  dev.off()
  
  legend_img <- magick::image_read(legend_file)
  legend_img <- magick::image_border(legend_img, color = "white", geometry = "20x20")
  
  legend_info <- magick::image_info(legend_img)
  legend_height <- legend_info$height
  
  top_margin <- 70
  side_margin <- 100
  gap_between_atlas_and_legend <- 70
  bottom_margin <- 70
  
  available_width <- final_width - (2 * side_margin)
  available_height <- final_height - top_margin - legend_height - gap_between_atlas_and_legend - bottom_margin
  
  atlas_img <- magick::image_resize(
    atlas_img,
    geometry = paste0(available_width, "x", available_height)
  )
  
  canvas <- magick::image_blank(
    width = final_width,
    height = final_height,
    color = "white"
  )
  
  canvas <- magick::image_composite(
    image = canvas,
    composite_image = atlas_img,
    gravity = "north",
    offset = paste0("+0+", top_margin)
  )
  
  canvas <- magick::image_composite(
    image = canvas,
    composite_image = legend_img,
    gravity = "south",
    offset = paste0("+0+", bottom_margin)
  )
  
  canvas <- magick::image_background(canvas, color = "white", flatten = TRUE)
  
  canvas <- magick::image_extent(
    canvas,
    geometry = paste0(final_width, "x", final_height),
    gravity = "center",
    color = "white"
  )
  
  if (tolower(format) %in% c("jpg", "jpeg")) {
    magick::image_write(canvas, path = file, format = "jpeg", quality = 95)
  } else {
    magick::image_write(canvas, path = file, format = "png")
  }
}

# ============================================================
# 6. Theme
# ============================================================

app_theme <- bs_theme(
  version = 5,
  bootswatch = "flatly",
  primary = "#1B7837",
  base_font = font_google("Source Sans 3"),
  heading_font = font_google("Source Sans 3")
)

# ============================================================
# 7. UI
# ============================================================

ui <- page_navbar(
  
  title = "Chickpea Gene Expression Atlas",
  theme = app_theme,
  bg = "#1B7837",
  inverse = TRUE,
  
  header = tags$head(
    tags$style(HTML("
      body {
        background-color: #f7f9fb;
      }
      .card {
        box-shadow: 0 2px 10px rgba(0,0,0,0.08);
        border-radius: 14px;
      }
      .atlas-title {
        font-size: 26px;
        font-weight: 700;
        color: #1B7837;
      }
      .small-note {
        color: #555;
        font-size: 14px;
      }
      .svg-box {
        background: white;
        border-radius: 14px;
        padding: 18px;
        border: 1px solid #e5e5e5;
        overflow-x: auto;
      }
      .download-btn {
        margin-bottom: 8px;
        width: 100%;
      }
      .about-section {
        line-height: 1.65;
      }
    "))
  ),
  
  nav_panel(
    "Expression Atlas",
    
    layout_sidebar(
      
      sidebar = sidebar(
        width = 340,
        
        h4("Gene search"),
        
        textInput(
          inputId = "single_gene_text",
          label = "Type gene ID",
          value = gene_choices[1],
          placeholder = "Example: Ca_v2.0_07355"
        ),
        
        tags$datalist(
          id = "gene_list",
          lapply(gene_choices, function(g) tags$option(value = g))
        ),
        
        tags$script(HTML("
          $(document).on('shiny:connected', function() {
            $('#single_gene_text').attr('list', 'gene_list');
          });
        ")),
        
        actionButton(
          inputId = "load_gene_atlas",
          label = "Load gene",
          class = "btn-success"
        ),
        
        br(), br(),
        
        uiOutput("atlas_gene_status"),
        
        sliderInput(
          inputId = "max_fpkm",
          label = "Maximum FPKM for atlas color scale",
          min = 1,
          max = 500,
          value = 100,
          step = 1
        ),
        
        selectInput(
          inputId = "palette",
          label = "Color palette",
          choices = names(palette_choices),
          selected = "Yellow-orange-red"
        ),
        
        selectInput(
          inputId = "stroke_mode",
          label = "Atlas outline/stroke",
          choices = c("No stroke", "Black stroke", "Gray stroke", "White stroke", "Custom color"),
          selected = "No stroke"
        ),
        
        conditionalPanel(
          condition = "input.stroke_mode == 'Custom color'",
          textInput(
            inputId = "custom_stroke_color",
            label = "Custom stroke color",
            value = "#333333",
            placeholder = "Example: #333333"
          )
        ),
        
        sliderInput(
          inputId = "stroke_width",
          label = "Stroke width",
          min = 0,
          max = 5,
          value = 1,
          step = 0.25
        ),
        
        hr(),
        
        h4("Download atlas"),
        
        downloadButton("download_atlas_svg", "Download atlas SVG", class = "download-btn"),
        downloadButton("download_atlas_png", "Download atlas PNG", class = "download-btn"),
        downloadButton("download_atlas_jpeg", "Download atlas JPEG", class = "download-btn")
      ),
      
      card(
        card_header(span(class = "atlas-title", "Expression eFP View")),
        p(class = "small-note", "Type a gene ID and click Load gene."),
        div(class = "svg-box", uiOutput("svg_output"))
      )
    )
  ),
  
  nav_panel(
    "Bar Graph",
    
    layout_sidebar(
      
      sidebar = sidebar(
        width = 330,
        
        h4("Gene and tissue options"),
        
        selectizeInput(
          inputId = "bar_gene",
          label = "Select gene",
          choices = gene_choices,
          selected = gene_choices[1],
          options = list(maxOptions = 5000)
        ),
        
        selectInput(
          inputId = "bar_tissue_group",
          label = "Tissue group",
          choices = names(tissue_groups),
          selected = "All tissues"
        ),
        
        checkboxInput(
          inputId = "bar_log",
          label = "Use log2(FPKM + 1)",
          value = TRUE
        ),
        
        hr(),
        
        h4("Download bar graph"),
        downloadButton("download_bar_png", "Download PNG", class = "download-btn"),
        downloadButton("download_bar_jpeg", "Download JPEG", class = "download-btn"),
        downloadButton("download_bar_pdf", "Download PDF", class = "download-btn")
      ),
      
      card(
        card_header("Expression across selected tissues"),
        plotlyOutput("barplot", height = "620px")
      )
    )
  ),
  
  nav_panel(
    "Mean / Median",
    
    layout_sidebar(
      
      sidebar = sidebar(
        width = 330,
        
        h4("Multiple gene selection"),
        
        selectizeInput(
          inputId = "summary_genes",
          label = "Select one or more genes",
          choices = gene_choices,
          selected = gene_choices[1],
          multiple = TRUE,
          options = list(maxOptions = 5000)
        ),
        
        textAreaInput(
          inputId = "custom_genes",
          label = "Or paste gene IDs here, one per line",
          placeholder = "Ca_v2.0_07355\nCa_v2.0_07371\nCa_v2.0_08859",
          rows = 6
        ),
        
        selectInput(
          inputId = "summary_tissue_group",
          label = "Tissue group",
          choices = names(tissue_groups),
          selected = "All tissues"
        ),
        
        checkboxInput(
          inputId = "summary_log",
          label = "Use log2(FPKM + 1)",
          value = TRUE
        ),
        
        hr(),
        
        h4("Download summary graph"),
        downloadButton("download_summary_png", "Download PNG", class = "download-btn"),
        downloadButton("download_summary_jpeg", "Download JPEG", class = "download-btn"),
        downloadButton("download_summary_csv", "Download summary CSV", class = "download-btn")
      ),
      
      card(
        card_header("Mean and median expression across selected genes"),
        plotlyOutput("summary_plot", height = "620px")
      )
    )
  ),
  
  nav_panel(
    "Heatmap",
    
    layout_sidebar(
      
      sidebar = sidebar(
        width = 330,
        
        h4("Heatmap gene options"),
        
        selectizeInput(
          inputId = "heatmap_genes",
          label = "Select genes",
          choices = gene_choices,
          selected = gene_choices[1:min(10, length(gene_choices))],
          multiple = TRUE,
          options = list(maxOptions = 5000)
        ),
        
        textAreaInput(
          inputId = "heatmap_custom_genes",
          label = "Or paste gene IDs here, one per line",
          placeholder = "Ca_v2.0_07355\nCa_v2.0_07371\nCa_v2.0_08859",
          rows = 6
        ),
        
        selectInput(
          inputId = "heatmap_tissue_group",
          label = "Tissue group",
          choices = names(tissue_groups),
          selected = "All tissues"
        ),
        
        selectInput(
          inputId = "heatmap_scale",
          label = "Heatmap value",
          choices = c("log2(FPKM + 1)", "Raw FPKM", "Row-scaled Z-score"),
          selected = "log2(FPKM + 1)"
        ),
        
        hr(),
        
        h4("Download heatmap"),
        downloadButton("download_heatmap_png", "Download PNG", class = "download-btn"),
        downloadButton("download_heatmap_jpeg", "Download JPEG", class = "download-btn"),
        downloadButton("download_heatmap_pdf", "Download PDF", class = "download-btn"),
        downloadButton("download_heatmap_csv", "Download heatmap CSV", class = "download-btn")
      ),
      
      card(
        card_header("Gene expression heatmap"),
        plotlyOutput("heatmap_plot", height = "720px")
      )
    )
  ),
  
  nav_panel(
    "GWAS Candidate Genes",
    
    layout_sidebar(
      
      sidebar = sidebar(
        width = 340,
        
        h4("SNP region search"),
        
        textInput(
          inputId = "gwas_snp_id",
          label = "SNP ID / marker name optional",
          value = "",
          placeholder = "Example: SNC_021162.2_67821180"
        ),
        
        textInput(
          inputId = "gwas_chr",
          label = "Chromosome",
          value = "3",
          placeholder = "Example: 3"
        ),
        
        numericInput(
          inputId = "gwas_snp_pos",
          label = "SNP position bp",
          value = 67821180,
          min = 1
        ),
        
        numericInput(
          inputId = "gwas_upstream",
          label = "Upstream distance bp",
          value = 50000,
          min = 0
        ),
        
        numericInput(
          inputId = "gwas_downstream",
          label = "Downstream distance bp",
          value = 50000,
          min = 0
        ),
        
        actionButton(
          inputId = "search_gwas_region",
          label = "Search candidate genes",
          class = "btn-success"
        ),
        
        hr(),
        
        downloadButton("download_gwas_candidates_csv", "Download CSV", class = "download-btn"),
        downloadButton("download_gwas_candidates_excel", "Download Excel", class = "download-btn")
      ),
      
      card(
        card_header("Candidate genes within selected SNP region"),
        uiOutput("gwas_region_text"),
        br(),
        DTOutput("gwas_candidate_table")
      )
    )
  ),
  
  nav_panel(
    "Expression Table",
    
    layout_sidebar(
      
      sidebar = sidebar(
        width = 330,
        
        h4("Table options"),
        
        selectizeInput(
          inputId = "table_genes",
          label = "Select genes",
          choices = gene_choices,
          selected = gene_choices[1],
          multiple = TRUE,
          options = list(maxOptions = 5000)
        ),
        
        textAreaInput(
          inputId = "table_custom_genes",
          label = "Or paste gene IDs here, one per line",
          rows = 6
        ),
        
        hr(),
        
        downloadButton("download_table_csv", "Download expression table CSV", class = "download-btn")
      ),
      
      card(
        card_header("Expression data table"),
        DTOutput("expr_table")
      )
    )
  ),
  
  nav_panel(
    "About",
    
    card(
      card_body(
        class = "about-section",
        
        div(
          style = "text-align:center; padding: 20px;",
          tags$img(
            src = "atlas_logo.png",
            style = "max-width: 220px; margin-bottom: 15px;"
          ),
          h2(
            "Chickpea Gene Expression Atlas",
            style = "color:#1B7837; font-weight:700; margin-bottom:5px;"
          ),
          h4(
            "Tennessee State University",
            style = "color:#333333; font-weight:600; margin-bottom:5px;"
          ),
          h5(
            "Version 1.0",
            style = "color:#666666; font-weight:500;"
          )
        ),
        
        hr(),
        
        h3("About the App", style = "color:#1B7837; font-weight:700;"),
        p("The Chickpea Gene Expression Atlas is an interactive web-based tool developed to visualize gene expression patterns across 32 chickpea tissues and developmental stages. The app allows users to search individual genes, view tissue-specific expression using an eFP-style atlas, generate expression bar plots, create heatmaps for selected genes, summarize mean and median expression, and identify candidate genes near GWAS-associated SNPs."),
        p("This resource is designed to support candidate gene prioritization, transcriptomic interpretation, and seed trait-related functional genomics research in chickpea."),
        
        hr(),
        
        h3("Authors", style = "color:#1B7837; font-weight:700;"),
        h4("Shubh Pravat Singh Yadav", style = "font-weight:700; margin-bottom:3px;"),
        p("College of Agriculture, Tennessee State University, Nashville, TN, USA"),
        p(HTML('Email: <a href="mailto:sushantpy8500@gmail.com">sushantpy8500@gmail.com</a><br>ORCID: <a href="https://orcid.org/0000-0003-3987-5616" target="_blank">0000-0003-3987-5616</a>')),
        
        h4("Kuber Shivashakarappa", style = "font-weight:700; margin-bottom:3px;"),
        p("College of Agriculture, Tennessee State University, Nashville, TN, USA"),
        p(HTML('Email: <a href="mailto:kshivash@my.tnstate.edu">kshivash@my.tnstate.edu</a>')),
        
        hr(),
        
        h3("Advisors", style = "color:#1B7837; font-weight:700;"),
        h4("Dr. Lyle Wallace", style = "font-weight:700; margin-bottom:3px;"),
        p("College of Agriculture, Tennessee State University, Nashville, TN, USA"),
        p(HTML('Email: <a href="mailto:lwalla10@tnstate.edu">lwalla10@tnstate.edu</a><br>ORCID: <a href="https://orcid.org/0000-0002-6985-5854" target="_blank">0000-0002-6985-5854</a>')),
        
        h4("Dr. Ali Taheri", style = "font-weight:700; margin-bottom:3px;"),
        p("College of Agriculture, Tennessee State University, Nashville, TN, USA"),
        p(HTML('Email: <a href="mailto:ali.taheri@tnstate.edu">ali.taheri@tnstate.edu</a><br>ORCID: <a href="https://orcid.org/0000-0002-4019-7982" target="_blank">0000-0002-4019-7982</a>')),
        
        hr(),
        
        h3("Affiliation", style = "color:#1B7837; font-weight:700;"),
        p(HTML('<b>College of Agriculture</b><br>Tennessee State University<br>3500 John A. Merritt Blvd<br>Nashville, TN 37209, USA')),
        
        hr(),
        
        h3("Credit and Data Source", style = "color:#1B7837; font-weight:700;"),
        p("The RNA-seq dataset used in this atlas was obtained from the publicly available chickpea transcriptome resource generated by Jain et al. (2022)."),
        p(HTML('<b>Reference:</b><br>Jain M, Bansal J, Rajkumar MS, Garg R. An integrated transcriptome mapping the regulatory network of coding and long non-coding RNAs provides a genomics resource in chickpea. <i>Communications Biology</i>. 2022;5:1106. PMID: 36261617.')),
        
        hr(),
        
        h3("Transcriptomic Analysis", style = "color:#1B7837; font-weight:700;"),
        p("Publicly available RNA-seq data from BioProject PRJNA622231 were used to examine expression profiles of candidate genes associated with chickpea seed traits. The dataset includes 32 tissues and developmental stages of chickpea cultivar ICC4958, including embryo, endosperm, seed coat, and seed samples collected at 5, 10, 20, and 30 days after pollination."),
        p("Transcript abundance was quantified using the CDC Frontier reference transcriptome. Gene-level expression values were summarized as FPKM and transformed as log2(FPKM + 1) for visualization across tissues."),
        
        hr(),
        
        h3("Website Development", style = "color:#1B7837; font-weight:700;"),
        p("The atlas website was developed as an R Shiny application using a custom SVG-based chickpea tissue diagram and gene-level expression matrix. Tissue-specific SVG elements were linked to expression values, allowing dynamic coloring of tissues based on selected gene expression."),
        p("The application also includes modules for bar plots, heatmaps, mean and median expression summaries, downloadable figures, and GWAS candidate gene searches using chromosome, SNP position, and upstream/downstream genomic windows."),
        
        hr(),
        
        div(
          style = "text-align:center; color:#666666; font-size:14px; padding-top:10px;",
          HTML("© Chickpea Gene Expression Atlas, Version 1.0<br>Tennessee State University")
        )
      )
    )
  )
)

# ============================================================
# 8. Server
# ============================================================

server <- function(input, output, session) {
  
  selected_palette <- reactive({
    palette_choices[[input$palette]]
  })
  
  selected_stroke_color <- reactive({
    if (input$stroke_mode == "Black stroke") {
      "#000000"
    } else if (input$stroke_mode == "Gray stroke") {
      "#333333"
    } else if (input$stroke_mode == "White stroke") {
      "#FFFFFF"
    } else if (input$stroke_mode == "Custom color") {
      input$custom_stroke_color
    } else {
      "none"
    }
  })
  
  current_atlas_gene <- reactiveVal(gene_choices[1])
  
  observeEvent(input$load_gene_atlas, {
    typed_gene <- trimws(input$single_gene_text)
    
    if (typed_gene %in% gene_choices) {
      current_atlas_gene(typed_gene)
    }
  })
  
  output$atlas_gene_status <- renderUI({
    typed_gene <- trimws(input$single_gene_text)
    
    if (typed_gene %in% gene_choices) {
      HTML("<span style='color:#1B7837; font-weight:bold;'>Gene found. Click Load gene.</span>")
    } else {
      HTML("<span style='color:#B2182B; font-weight:bold;'>Gene not found in expression matrix.</span>")
    }
  })
  
  atlas_svg_text <- reactive({
    
    gene_row <- expr %>%
      filter(.data[[gene_col]] == current_atlas_gene())
    
    if (nrow(gene_row) == 0) {
      return("<b>Gene not found.</b>")
    }
    
    color_svg(
      gene_row = gene_row,
      max_fpkm = input$max_fpkm,
      pal_colors = selected_palette(),
      stroke_mode = input$stroke_mode,
      stroke_color = selected_stroke_color(),
      stroke_width = input$stroke_width
    )
  })
  
  output$svg_output <- renderUI({
    HTML(atlas_svg_text())
  })
  
  output$download_atlas_svg <- downloadHandler(
    filename = function() {
      paste0(current_atlas_gene(), "_expression_atlas.svg")
    },
    content = function(file) {
      writeLines(atlas_svg_text(), file)
    }
  )
  
  output$download_atlas_png <- downloadHandler(
    filename = function() {
      paste0(current_atlas_gene(), "_expression_atlas.png")
    },
    content = function(file) {
      save_atlas_with_legend(
        svg_text = atlas_svg_text(),
        file = file,
        format = "png",
        final_width = 4500,
        final_height = 3250,
        max_fpkm = input$max_fpkm,
        pal_colors = selected_palette()
      )
    }
  )
  
  output$download_atlas_jpeg <- downloadHandler(
    filename = function() {
      paste0(current_atlas_gene(), "_expression_atlas.jpeg")
    },
    content = function(file) {
      save_atlas_with_legend(
        svg_text = atlas_svg_text(),
        file = file,
        format = "jpeg",
        final_width = 4500,
        final_height = 3250,
        max_fpkm = input$max_fpkm,
        pal_colors = selected_palette()
      )
    }
  )
  
  # ----------------------------------------------------------
  # Bar plot
  # ----------------------------------------------------------
  
  bar_data <- reactive({
    selected_tissues <- tissue_groups[[input$bar_tissue_group]]
    
    expr %>%
      filter(.data[[gene_col]] == input$bar_gene) %>%
      select(all_of(gene_col), all_of(selected_tissues)) %>%
      pivot_longer(
        cols = all_of(selected_tissues),
        names_to = "Tissue",
        values_to = "FPKM"
      ) %>%
      mutate(
        Tissue = factor(Tissue, levels = selected_tissues),
        FPKM = as.numeric(FPKM),
        Log2_FPKM = log2(FPKM + 1)
      )
  })
  
  make_bar_plot <- function() {
    df <- bar_data()
    
    yvar <- ifelse(input$bar_log, "Log2_FPKM", "FPKM")
    ylab <- ifelse(input$bar_log, "log2(FPKM + 1)", "FPKM")
    
    ggplot(df, aes(x = reorder(Tissue, .data[[yvar]]), y = .data[[yvar]])) +
      geom_col(fill = "#2C7FB8", width = 0.75) +
      coord_flip() +
      theme_bw(base_size = 14) +
      theme(
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 11),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      ) +
      labs(
        title = paste("Expression profile of", input$bar_gene),
        x = "Tissue",
        y = ylab
      )
  }
  
  output$barplot <- renderPlotly({
    df <- bar_data()
    yvar <- ifelse(input$bar_log, "Log2_FPKM", "FPKM")
    ylab <- ifelse(input$bar_log, "log2(FPKM + 1)", "FPKM")
    
    p <- ggplot(
      df,
      aes(
        x = reorder(Tissue, .data[[yvar]]),
        y = .data[[yvar]],
        text = paste(
          "Tissue:", Tissue,
          "<br>FPKM:", round(FPKM, 3),
          "<br>", ylab, ":", round(.data[[yvar]], 3)
        )
      )
    ) +
      geom_col(fill = "#2C7FB8", width = 0.75) +
      coord_flip() +
      theme_bw(base_size = 14) +
      theme(
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", size = 16),
        axis.text.y = element_text(size = 11)
      ) +
      labs(
        title = paste("Expression profile of", input$bar_gene),
        x = "Tissue",
        y = ylab
      )
    
    ggplotly(p, tooltip = "text")
  })
  
  output$download_bar_png <- downloadHandler(
    filename = function() paste0(input$bar_gene, "_barplot.png"),
    content = function(file) ggsave(file, make_bar_plot(), width = 10, height = 8, dpi = 600, bg = "white")
  )
  
  output$download_bar_jpeg <- downloadHandler(
    filename = function() paste0(input$bar_gene, "_barplot.jpeg"),
    content = function(file) ggsave(file, make_bar_plot(), width = 10, height = 8, dpi = 600, bg = "white")
  )
  
  output$download_bar_pdf <- downloadHandler(
    filename = function() paste0(input$bar_gene, "_barplot.pdf"),
    content = function(file) ggsave(file, make_bar_plot(), width = 10, height = 8, bg = "white")
  )
  
  # ----------------------------------------------------------
  # Mean / Median
  # ----------------------------------------------------------
  
  summary_data <- reactive({
    genes <- parse_gene_input(input$summary_genes, input$custom_genes)
    selected_tissues <- tissue_groups[[input$summary_tissue_group]]
    
    validate(need(length(genes) > 0, "Please select or paste at least one valid gene ID."))
    
    expr %>%
      filter(.data[[gene_col]] %in% genes) %>%
      select(all_of(gene_col), all_of(selected_tissues)) %>%
      pivot_longer(
        cols = all_of(selected_tissues),
        names_to = "Tissue",
        values_to = "FPKM"
      ) %>%
      mutate(
        FPKM = as.numeric(FPKM),
        Value = ifelse(input$summary_log, log2(FPKM + 1), FPKM)
      ) %>%
      group_by(Tissue) %>%
      summarise(
        Mean = mean(Value, na.rm = TRUE),
        Median = median(Value, na.rm = TRUE),
        Gene_count = n_distinct(.data[[gene_col]]),
        .groups = "drop"
      ) %>%
      mutate(Tissue = factor(Tissue, levels = selected_tissues))
  })
  
  make_summary_plot <- function() {
    df <- summary_data() %>%
      pivot_longer(cols = c("Mean", "Median"), names_to = "Statistic", values_to = "Expression")
    
    ylab <- ifelse(input$summary_log, "log2(FPKM + 1)", "FPKM")
    
    ggplot(df, aes(x = Tissue, y = Expression, fill = Statistic)) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      coord_flip() +
      theme_bw(base_size = 14) +
      theme(
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", size = 16),
        legend.position = "top",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      ) +
      labs(
        title = "Mean and median expression across selected genes",
        x = "Tissue",
        y = ylab,
        fill = "Statistic"
      )
  }
  
  output$summary_plot <- renderPlotly({
    df <- summary_data() %>%
      pivot_longer(cols = c("Mean", "Median"), names_to = "Statistic", values_to = "Expression")
    
    ylab <- ifelse(input$summary_log, "log2(FPKM + 1)", "FPKM")
    
    p <- ggplot(
      df,
      aes(
        x = Tissue,
        y = Expression,
        fill = Statistic,
        text = paste(
          "Tissue:", Tissue,
          "<br>Statistic:", Statistic,
          "<br>Expression:", round(Expression, 3),
          "<br>Genes:", Gene_count
        )
      )
    ) +
      geom_col(position = position_dodge(width = 0.8), width = 0.7) +
      coord_flip() +
      theme_bw(base_size = 14) +
      theme(
        panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", size = 16),
        legend.position = "top"
      ) +
      labs(
        title = "Mean and median expression across selected genes",
        x = "Tissue",
        y = ylab,
        fill = "Statistic"
      )
    
    ggplotly(p, tooltip = "text")
  })
  
  output$download_summary_png <- downloadHandler(
    filename = function() "mean_median_expression.png",
    content = function(file) ggsave(file, make_summary_plot(), width = 11, height = 8, dpi = 600, bg = "white")
  )
  
  output$download_summary_jpeg <- downloadHandler(
    filename = function() "mean_median_expression.jpeg",
    content = function(file) ggsave(file, make_summary_plot(), width = 11, height = 8, dpi = 600, bg = "white")
  )
  
  output$download_summary_csv <- downloadHandler(
    filename = function() "mean_median_expression_summary.csv",
    content = function(file) write.csv(summary_data(), file, row.names = FALSE)
  )
  
  # ----------------------------------------------------------
  # Heatmap
  # ----------------------------------------------------------
  
  heatmap_data <- reactive({
    genes <- parse_gene_input(input$heatmap_genes, input$heatmap_custom_genes)
    selected_tissues <- tissue_groups[[input$heatmap_tissue_group]]
    
    validate(need(length(genes) > 0, "Please select or paste at least one valid gene ID."))
    
    df <- expr %>%
      filter(.data[[gene_col]] %in% genes) %>%
      select(all_of(gene_col), all_of(selected_tissues)) %>%
      pivot_longer(
        cols = all_of(selected_tissues),
        names_to = "Tissue",
        values_to = "FPKM"
      ) %>%
      mutate(
        FPKM = as.numeric(FPKM),
        Value = case_when(
          input$heatmap_scale == "Raw FPKM" ~ FPKM,
          input$heatmap_scale == "log2(FPKM + 1)" ~ log2(FPKM + 1),
          TRUE ~ log2(FPKM + 1)
        )
      )
    
    if (input$heatmap_scale == "Row-scaled Z-score") {
      df <- df %>%
        group_by(.data[[gene_col]]) %>%
        mutate(Value = as.numeric(scale(log2(FPKM + 1)))) %>%
        ungroup()
    }
    
    df %>%
      mutate(
        Tissue = factor(Tissue, levels = selected_tissues),
        Gene = factor(.data[[gene_col]], levels = rev(unique(genes)))
      )
  })
  
  make_heatmap_plot <- function() {
    df <- heatmap_data()
    
    ggplot(df, aes(x = Tissue, y = Gene, fill = Value)) +
      geom_tile(color = "white", linewidth = 0.25) +
      scale_fill_gradientn(
        colors = c("#f7fbff", "#ffffb2", "#fd8d3c", "#bd0026"),
        name = input$heatmap_scale
      ) +
      theme_bw(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 16),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA)
      ) +
      labs(
        title = "Gene expression heatmap",
        x = "Tissue",
        y = "Gene"
      )
  }
  
  output$heatmap_plot <- renderPlotly({
    df <- heatmap_data()
    
    p <- ggplot(
      df,
      aes(
        x = Tissue,
        y = Gene,
        fill = Value,
        text = paste(
          "Gene:", Gene,
          "<br>Tissue:", Tissue,
          "<br>FPKM:", round(FPKM, 3),
          "<br>Value:", round(Value, 3)
        )
      )
    ) +
      geom_tile(color = "white", linewidth = 0.25) +
      scale_fill_gradientn(
        colors = c("#f7fbff", "#ffffb2", "#fd8d3c", "#bd0026"),
        name = input$heatmap_scale
      ) +
      theme_bw(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 9),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 16)
      ) +
      labs(
        title = "Gene expression heatmap",
        x = "Tissue",
        y = "Gene"
      )
    
    ggplotly(p, tooltip = "text")
  })
  
  output$download_heatmap_png <- downloadHandler(
    filename = function() "gene_expression_heatmap.png",
    content = function(file) {
      n_genes <- length(unique(heatmap_data()[[gene_col]]))
      h <- max(6, min(20, n_genes * 0.35))
      ggsave(file, make_heatmap_plot(), width = 13, height = h, dpi = 600, bg = "white")
    }
  )
  
  output$download_heatmap_jpeg <- downloadHandler(
    filename = function() "gene_expression_heatmap.jpeg",
    content = function(file) {
      n_genes <- length(unique(heatmap_data()[[gene_col]]))
      h <- max(6, min(20, n_genes * 0.35))
      ggsave(file, make_heatmap_plot(), width = 13, height = h, dpi = 600, bg = "white")
    }
  )
  
  output$download_heatmap_pdf <- downloadHandler(
    filename = function() "gene_expression_heatmap.pdf",
    content = function(file) {
      n_genes <- length(unique(heatmap_data()[[gene_col]]))
      h <- max(6, min(20, n_genes * 0.35))
      ggsave(file, make_heatmap_plot(), width = 13, height = h, bg = "white")
    }
  )
  
  output$download_heatmap_csv <- downloadHandler(
    filename = function() "heatmap_expression_values.csv",
    content = function(file) write.csv(heatmap_data(), file, row.names = FALSE)
  )
  
  # ----------------------------------------------------------
  # GWAS candidate gene finder
  # ----------------------------------------------------------
  
  gwas_candidates <- eventReactive(input$search_gwas_region, {
    chr_input <- normalize_chr(input$gwas_chr)
    snp_pos <- as.numeric(input$gwas_snp_pos)
    upstream <- as.numeric(input$gwas_upstream)
    downstream <- as.numeric(input$gwas_downstream)
    
    region_start <- max(1, snp_pos - upstream)
    region_end <- snp_pos + downstream
    
    gene_annot %>%
      mutate(Chr_clean = normalize_chr(Chr)) %>%
      filter(
        Chr_clean == chr_input,
        End >= region_start,
        Start <= region_end
      ) %>%
      mutate(
        SNP_ID = ifelse(input$gwas_snp_id == "", NA, input$gwas_snp_id),
        SNP_position = snp_pos,
        Region_start = region_start,
        Region_end = region_end,
        Gene_midpoint = (Start + End) / 2,
        Distance_from_SNP_bp = round(abs(Gene_midpoint - snp_pos), 0),
        Direction = case_when(
          End < snp_pos ~ "Upstream",
          Start > snp_pos ~ "Downstream",
          TRUE ~ "SNP within/overlapping gene"
        ),
        Present_in_expression_matrix = ifelse(Gene_ID %in% expr[[gene_col]], "Yes", "No")
      ) %>%
      select(
        SNP_ID,
        Chr,
        SNP_position,
        Region_start,
        Region_end,
        Gene_ID,
        Gene_raw,
        Start,
        End,
        Strand,
        Gene_length_bp,
        Distance_from_SNP_bp,
        Direction,
        Annotation,
        Dbxref,
        Present_in_expression_matrix,
        everything()
      ) %>%
      arrange(Distance_from_SNP_bp)
  })
  
  output$gwas_region_text <- renderUI({
    req(input$search_gwas_region)
    
    chr_input <- normalize_chr(input$gwas_chr)
    snp_pos <- as.numeric(input$gwas_snp_pos)
    upstream <- as.numeric(input$gwas_upstream)
    downstream <- as.numeric(input$gwas_downstream)
    
    region_start <- max(1, snp_pos - upstream)
    region_end <- snp_pos + downstream
    
    n_genes <- nrow(gwas_candidates())
    
    HTML(paste0(
      "<b>Search region:</b> Chr", chr_input,
      ": ", format(region_start, big.mark = ","),
      " - ", format(region_end, big.mark = ","),
      "<br><b>SNP position:</b> ",
      format(snp_pos, big.mark = ","),
      "<br><b>Candidate genes found:</b> ",
      n_genes
    ))
  })
  
  output$gwas_candidate_table <- renderDT({
    datatable(
      gwas_candidates(),
      options = list(
        pageLength = 20,
        scrollX = TRUE
      )
    )
  })
  
  output$download_gwas_candidates_csv <- downloadHandler(
    filename = function() {
      paste0("GWAS_candidate_genes_Chr", input$gwas_chr, "_", input$gwas_snp_pos, ".csv")
    },
    content = function(file) {
      write.csv(gwas_candidates(), file, row.names = FALSE)
    }
  )
  
  output$download_gwas_candidates_excel <- downloadHandler(
    filename = function() {
      paste0("GWAS_candidate_genes_Chr", input$gwas_chr, "_", input$gwas_snp_pos, ".xlsx")
    },
    content = function(file) {
      openxlsx::write.xlsx(gwas_candidates(), file, overwrite = TRUE)
    }
  )
  
  # ----------------------------------------------------------
  # Expression table
  # ----------------------------------------------------------
  
  table_data <- reactive({
    genes <- parse_gene_input(input$table_genes, input$table_custom_genes)
    
    validate(need(length(genes) > 0, "Please select or paste at least one valid gene ID."))
    
    expr %>%
      filter(.data[[gene_col]] %in% genes)
  })
  
  output$expr_table <- renderDT({
    datatable(
      table_data(),
      options = list(
        pageLength = 20,
        scrollX = TRUE
      )
    )
  })
  
  output$download_table_csv <- downloadHandler(
    filename = function() "selected_gene_expression_table.csv",
    content = function(file) write.csv(table_data(), file, row.names = FALSE)
  )
}

# ============================================================
# 9. Run app
# ============================================================

shinyApp(ui, server)