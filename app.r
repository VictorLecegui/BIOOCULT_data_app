library(shiny)
library(dplyr)
library(tidyr)
library(leaflet)

# =========================
# 1. LOAD & PREPARE DATA
# =========================

metadata2 <- read.csv("data/Metadata_FW.csv")  |> 
  separate(Coords, into = c("lat", "lon"), sep = ",")  |> 
  mutate(
    lon = trimws(lon),
    lat = as.numeric(lat),
    lon = as.numeric(lon)
  )


db <- readRDS("data/OCS_normalized_2026.rds") |> 
  filter(!(obitag_rank %in% c("no rank", "domain", "kingdom", "cellular root"))) |> 
  filter(!(is.na(domain_name) | domain_name=="Bacteria")) |> 
  pivot_longer(cols = c(27:605), names_to = "Plot_ID", values_to = "N_reads") |> 
              filter(N_reads > 0)  |>  
  select(id, Plot_ID, N_reads, 
        phylum_name, class_name, order_name, 
        family_name, genus_name, species_name, 
        scientific_name) |> 
  rename(Phylum = phylum_name, 
          Class = class_name, 
          Order = order_name, 
          Family = family_name, 
          Genus = genus_name, 
          Species = species_name, 
          Taxon = scientific_name) |> 
         left_join(metadata2, by = "Plot_ID") %>%
  mutate(across(Phylum:Taxon, ~replace_na(., "Unknown"))) %>%
  distinct(id, Plot_ID, .keep_all = TRUE)


# Taxonomy tree
taxo_tree <- db %>%
  distinct(Phylum, Class, Order, Family, Genus, Species)

# =========================
# 2. UI
# =========================
ui <- fluidPage(
  titlePanel("MOTU Explorer"),

  sidebarLayout(
    sidebarPanel(
      selectizeInput("phylum", "Phylum",
                     choices = sort(unique(db$Phylum)),
                     multiple = TRUE,
                     options = list(create = TRUE, placeholder = "Type or select")),

      selectizeInput("class", "Class", 
                    choices = sort(unique(db$Class)), multiple = TRUE,
                     options = list(create = TRUE, placeholder = "Type or select")),

      selectizeInput("order", "Order", 
                        choices = sort(unique(db$Order)) , multiple = TRUE,
                     options = list(create = TRUE, placeholder = "Type or select")),

      selectizeInput("family", "Family", 
                    choices = sort(unique(db$Family)), multiple = TRUE,
                     options = list(create = TRUE, placeholder = "Type or select")),

      selectizeInput("genus", "Genus", 
                    choices = sort(unique(db$Genus)) , multiple = TRUE,
                     options = list(create = TRUE, placeholder = "Type or select")),

      selectizeInput("species", "Species", 
                    choices = sort(unique(db$Species)), multiple = TRUE,
                     options = list(create = TRUE, placeholder = "Type or select")),
      selectizeInput("taxon", "Taxon", 
                    choices = sort(unique(db$Taxon)), multiple = TRUE,
                     options = list(create = TRUE, placeholder = "Type or select")),

      hr(),

      downloadButton("download_data", "Download filtered plots")
    ),

    mainPanel(
      leafletOutput("map", height = 500),
      br(),
      tableOutput("table")
    )
  )
)
# =========================
# 3. SERVER
# =========================

server <- function(input, output, session) {

  # ---- Helper: flexible matching ----
  matches_input <- function(x, input_vec) {
    if (is.null(input_vec) || length(input_vec) == 0) return(rep(TRUE, length(x)))
    grepl(paste(input_vec, collapse = "|"), x, ignore.case = TRUE)
  }

  # ---- Cascading filters ----

  observeEvent(input$phylum, {
    filtered <- taxo_tree %>%
      filter(matches_input(Phylum, input$phylum))

    updateSelectizeInput(session, "class",
                         choices = sort(unique(filtered$Class)),
                         selected = character(0))
  })

  observeEvent(input$class, {
    filtered <- taxo_tree %>%
      filter(
        matches_input(Phylum, input$phylum),
        matches_input(Class, input$class)
      )

    updateSelectizeInput(session, "order",
                         choices = sort(unique(filtered$Order)),
                         selected = character(0))
  })

  observeEvent(input$order, {
    filtered <- taxo_tree %>%
      filter(
        matches_input(Phylum, input$phylum),
        matches_input(Class, input$class),
        matches_input(Order, input$order)
      )

    updateSelectizeInput(session, "family",
                         choices = sort(unique(filtered$Family)),
                         selected = character(0))
  })

  observeEvent(input$family, {
    filtered <- taxo_tree %>%
      filter(
        matches_input(Phylum, input$phylum),
        matches_input(Class, input$class),
        matches_input(Order, input$order),
        matches_input(Family, input$family)
      )

    updateSelectizeInput(session, "genus",
                         choices = sort(unique(filtered$Genus)),
                         selected = character(0))
  })

  observeEvent(input$genus, {
    filtered <- taxo_tree %>%
      filter(
        matches_input(Phylum, input$phylum),
        matches_input(Class, input$class),
        matches_input(Order, input$order),
        matches_input(Family, input$family),
        matches_input(Genus, input$genus)
      )

    updateSelectizeInput(session, "species",
                         choices = sort(unique(filtered$Species)),
                         selected = character(0))
  })

  # ---- Final filtering ----

  filtered_db <- reactive({
    db %>%
      filter(
        matches_input(Phylum, input$phylum),
        matches_input(Class, input$class),
        matches_input(Order, input$order),
        matches_input(Family, input$family),
        matches_input(Genus, input$genus),
        matches_input(Species, input$species), 
        matches_input(Taxon, input$taxon)
      )
  })

  # ---- Extract plots ----

  plots <- reactive({
    filtered_db() %>%
      distinct(Plot_ID, LOC, lat, lon, Elev, pH, CE, WTD_mean,
               Veg_height_mean, Pertorb, Tmit_anual, Pluvio_anual,
               pv_groups, bry_groups, PV_N_sp, BRY_N_sp)
  })

  # ---- Map with basemap switch ----

  output$map <- renderLeaflet({
    leaflet(plots()) %>%

      addProviderTiles("OpenStreetMap", group = "OSM") %>%
      addProviderTiles("Esri.WorldImagery", group = "Orto") %>%

      addCircleMarkers(
        ~lon, ~lat,
        popup = ~paste0("<b>", Plot_ID, "</b><br>",
                        "Localitat: ", LOC, "<br>",
                        "Elev: ", Elev),
        group = "Plots"
      ) %>%

      addLayersControl(
        baseGroups = c("OSM", "Orto"),
        overlayGroups = c("Plots"),
        options = layersControlOptions(collapsed = FALSE)
      )
  })

  # ---- Table ----

  output$table <- renderTable({
    plots()
  })

  # ---- Download ----

  output$download_data <- downloadHandler(
    filename = function() {
      "filtered_plots.csv"
    },
    content = function(file) {
      write.csv(filtered_db(), file, row.names = FALSE)
    }
  )
}

# =========================
# 4. RUN APP
# =========================

shinyApp(ui, server)
