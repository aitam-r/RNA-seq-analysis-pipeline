# Sources ----
source("ui_gsu.R", local = TRUE)
source("ui_deseq2.R", local = TRUE)
source("ui_wgcna.R", local = TRUE)

# UI -----------------------------------------------------------------------------

ui <- dashboardPage(
  dashboardHeader(title = "RNAinsights"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("General Set-up", tabName = "gsu"),
      menuItem("DESeq2", tabName = "deseq2",
               menuSubItem("DESeq2 Set-up", tabName = "deseq2_su"),
               menuSubItem("Exploration", tabName = "explo"),
               menuSubItem("Many plots", tabName = "plots"),
               menuSubItem("DEG table", tabName = "deg"),
               menuSubItem("Single gene counts", tabName = "sgc")),
      menuItem("Co-expression Network", tabName = "wgcna",
               menuSubItem("WGCNA Set-up", tabName = "wgcna_su"),
               menuSubItem("Modules' merging", tabName = "mod_merge"),
               menuSubItem("Modules' counts", tabName = "mod_count"),
               menuSubItem("Enrichment", tabName = "enr"))
    )
  ),
  dashboardBody(
    shinyDashboardThemes(theme = "poor_mans_flatly"),
    
    tabItems(
      tab_gsu,
      tab_deseq2_su,
      tab_explo,
      tab_plots,
      tab_deg,
      tab_sgc,
      tab_wgcna_su,
      tab_mod_merge,
      tab_mod_count,
      tab_enr
    )
  )
)

