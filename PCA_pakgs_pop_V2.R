## Relevance of selected packages
## based on downloads from CRAN

## resources for download stats
## -- https://stackoverflow.com/questions/75411383/how-do-i-know-if-an-r-package-is-commonly-used
## -- https://bioconductor.org/packages/stats/
## existing applets
## -- https://dgrtwo.shinyapps.io/cranview/

## Notice that Gower's package is unavailable on CRAN and must be retrieved from the book's website: https://www.wiley.com/legacy/wileychi/gower/material.html

## Own function to reshape bioconductor download data ####

changeBiocondate <- function(pkgName, pkgsData, fromDate, toDate, chosen_colnames){
  ## step1: remove "all months"
  idx_allmonths <- which(pkgsData$Month == "all")
  downloads_pkg <- pkgsData[-idx_allmonths,]
  ## change dates --  chosen
  ## -- see https://stackoverflow.com/questions/39762543/convert-character-month-name-to-date-time-object
  month_to_number <- match(downloads_pkg$Month, month.abb)
  new_date <- paste("1", month_to_number, downloads_pkg$Year, sep="/")
  downloads_pkg$date <- as.POSIXct(new_date, format = "%d/%m/%Y", tz = "GMT")
  to_keep <- c(which(downloads_pkg$date >= fromDate & downloads_pkg$date <= toDate))
  
  downloads_pkg <- downloads_pkg[to_keep, c("date", "Nb_of_downloads")]
  downloads_pkg$package <- pkgName
  colnames(downloads_pkg) <- chosen_colnames
  return(downloads_pkg)
}


## Load data and packages                             ####
library(cranlogs)
mypath <- "./"
Out_path <- paste0(mypath, "Analysis_OUT")   

## get my extended review from spreadsheet (extend my initial draft after review)
all_other_packs <- read.csv(paste0(mypath, "packages_from_search.csv"), header = T)
all_other_packs_cran <- all_other_packs[which(all_other_packs$Repository == "CRAN"),]
all_other_packs_bioc <- all_other_packs[which(all_other_packs$Repository == "Bioconductor"),]

##                                                    ####

## Part 1: Get donwload stats                         ####
## -- Get monthly downloads, cran pkgs                ####
all_CRAN_downloads <- cran_downloads(packages = all_other_packs_cran$package, 
                                      from = "2023-03-01",
                                      to = "2025-03-01")
chosen_colnames <- colnames(all_CRAN_downloads)                                   # will be useful later




## -- Get monthly downloads, bioconductor pkgs        ####

## Download stats from url                           
bioc_stats_1 <- url("https://www.bioconductor.org/packages/stats/bioc/pcaMethods/pcaMethods_stats.tab")
bioc_stats_2 <- url("https://www.bioconductor.org/packages/stats/bioc/PCAtools/PCAtools_stats.tab")
bioc_stats_3 <- url("https://www.bioconductor.org/packages/stats/bioc/scPCA/scPCA_stats.tab")
bioc_stats_4 <- url("https://www.bioconductor.org/packages/stats/bioc/regionalpcs/regionalpcs_stats.tab")
bioc_stats_5 <- url("https://www.bioconductor.org/packages/stats/bioc/SNPRelate/SNPRelate_stats.tab")
bioc_stats_6 <- url("https://www.bioconductor.org/packages/stats/bioc/ropls/ropls_stats.tab")
bioc_stats_7 <- url("https://www.bioconductor.org/packages/stats/bioc/pcaExplorer/pcaExplorer_stats.tab")

downloads_pcaMethods = read.table(bioc_stats_1, header=TRUE)             # https://stackoverflow.com/questions/49155440/downloading-tab-delimited-file-from-internet
downloads_pcaTools = read.table(bioc_stats_2,header=TRUE)
downloads_scPCA = read.table(bioc_stats_3, header=TRUE) 
downloads_regionalpcs = read.table(bioc_stats_4, header=TRUE) 
downloads_SNPRelate = read.table(bioc_stats_5, header=TRUE) 
downloads_ropls = read.table(bioc_stats_6, header=TRUE) 
downloads_pcaExplorer = read.table(bioc_stats_7, header=TRUE) 

pkgData_biocon <- list(downloads_pcaMethods,
                       downloads_pcaTools,
                       downloads_scPCA,
                       downloads_regionalpcs,
                       downloads_SNPRelate,
                       downloads_ropls,
                       downloads_pcaExplorer
)

## need to name in the order the files have been downloaded
pkgName_biocon <- c("pcaMethods", 
                    "PCAtools",
                    "scPCA",
                    "regionalpcs",
                    "SNPRelate",
                    "ropls",
                    "pcaExplorer")



## create list then call my own function to reshape bioconductor data's format
names(pkgData_biocon) <- pkgName_biocon 
fromDate <- "2023-03-01" 
toDate <- "2025-03-01"
reshaped_list_biocDownload <- lapply(pkgName_biocon, function(x){
  pkgName <- x
  pkgsData <- pkgData_biocon[[x]]
  temp_data <- changeBiocondate(pkgName = pkgName, 
                                pkgsData = pkgsData,
                                fromDate = fromDate,
                                toDate = toDate,
                                chosen_colnames = chosen_colnames)
  return(temp_data)
})
all_bioc_downloads <- do.call("rbind", reshaped_list_biocDownload )


## -- Compute total downloads per package             ####
all_bioc_downloads_total <- aggregate(all_bioc_downloads$count, by = list(all_bioc_downloads$package), FUN=sum)
all_CRAN_downloads_total <- aggregate(all_CRAN_downloads$count, by = list(all_CRAN_downloads$package), FUN=sum)

## update original dataset with counts
pkg_idx_bioc <- match(all_bioc_downloads_total$Group.1, all_other_packs_bioc$package)                 # find pkg name poisition in original database
pkgs_idx_CRAN <- match(all_CRAN_downloads_total$Group.1, all_other_packs_cran$package)
all_other_packs_bioc$download_count <- NA
all_other_packs_cran$download_count <- NA
all_other_packs_bioc$download_count[pkg_idx_bioc] <- all_bioc_downloads_total$x
all_other_packs_cran$download_count[pkgs_idx_CRAN] <- all_CRAN_downloads_total$x





## Part 2: rank pkgs                                  ####
## -- Consolidate and plot download stats
all_pkgs_total_downloads <- rbind(all_other_packs_cran, all_other_packs_bioc)

all_pkgs_total_downloads <- all_pkgs_total_downloads[order(all_pkgs_total_downloads$download_count, decreasing = T),]
all_pkgs_total_downloads_TOP20 <- all_pkgs_total_downloads[1:20,]

## top 10, PCA only, for subplot
all_pkgs_subset1_downloads <- all_pkgs_total_downloads[which(all_pkgs_total_downloads$Search_type == "PCA_search"),]
all_pkgs_subset1_downloads_TOP10 <- all_pkgs_subset1_downloads[1:10,]

## top 10, biplot only, for subplot
all_pkgs_subset2_downloads <- all_pkgs_total_downloads[which(all_pkgs_total_downloads$Search_type == "BIPLOT_search"),]
all_pkgs_subset2_downloads_TOP10 <- all_pkgs_subset2_downloads[1:10,]

## Note the counterintuitive sorting due to the weird behaviour of the horizontal barplot
par(mar = c(4,6,2,2))
barplot(height= all_pkgs_total_downloads_TOP20$download_count, 
        names= all_pkgs_total_downloads_TOP20$package,
        #col = col.colors[all_pkgs_total_downloads_TOP20$repo],
        horiz = T,
        las = 1)                                                               # labels rotation: https://r-graph-gallery.com/210-custom-barplot-layout.html; https://stackoverflow.com/questions/10286473/rotating-x-axis-labels-in-r-for-barplot


## Part 3: Plot download data for final subset selection                  ####



boxplot(count ~ package,
        data = all_pkgs_total_downloads[which(all_pkgs_total_downloads$Color_coding != "Excluded"),] ,
        xlab = "package name",
        ylab = "monthly downloads")




## Part 3: use  plotly tp visualise for paper ####

all_pkgs_total_downloads_TOP20$Color_coding2 <- NA
all_pkgs_subset1_downloads_TOP10$Color_coding2 <- NA 
all_pkgs_subset2_downloads_TOP10$Color_coding2 <- NA 

all_pkgs_total_downloads_TOP20$Color_coding2[which(all_pkgs_total_downloads_TOP20$Color_coding == "Generalist")] <- "incl.-Generalist"
all_pkgs_total_downloads_TOP20$Color_coding2[which(all_pkgs_total_downloads_TOP20$Color_coding == "PCA")] <- "incl.-specialist-PCA"
all_pkgs_total_downloads_TOP20$Color_coding2[which(all_pkgs_total_downloads_TOP20$Color_coding == "Biplot")] <- "incl.-specialist-Biplot"
all_pkgs_total_downloads_TOP20$Color_coding2[which(all_pkgs_total_downloads_TOP20$Color_coding == "Excluded")] <- "Excluded"

all_pkgs_subset1_downloads_TOP10$Color_coding2[which(all_pkgs_subset1_downloads_TOP10$Color_coding == "Generalist")] <- "incl.-Generalist"
all_pkgs_subset1_downloads_TOP10$Color_coding2[which(all_pkgs_subset1_downloads_TOP10$Color_coding == "PCA")] <- "incl.-specialist-PCA"
all_pkgs_subset1_downloads_TOP10$Color_coding2[which(all_pkgs_subset1_downloads_TOP10$Color_coding == "Biplot")] <- "incl.-specialist-Biplot"
all_pkgs_subset1_downloads_TOP10$Color_coding2[which(all_pkgs_subset1_downloads_TOP10$Color_coding == "Excluded")] <- "Excluded"

all_pkgs_subset2_downloads_TOP10$Color_coding2[which(all_pkgs_subset2_downloads_TOP10$Color_coding == "Generalist")] <- "incl.-Generalist"
all_pkgs_subset2_downloads_TOP10$Color_coding2[which(all_pkgs_subset2_downloads_TOP10$Color_coding == "PCA")] <- "incl.-specialist-PCA"
all_pkgs_subset2_downloads_TOP10$Color_coding2[which(all_pkgs_subset2_downloads_TOP10$Color_coding == "Biplot")] <- "incl.-specialist-Biplot"
all_pkgs_subset2_downloads_TOP10$Color_coding2[which(all_pkgs_subset2_downloads_TOP10$Color_coding == "Excluded")] <- "Excluded"




## thread: https://plotly.com/r/horizontal-bar-charts/
## for the defaul color names: https://community.plotly.com/t/plotly-colours-list/11730/2
## on combining plots: https://plotly.com/r/subplots/
## on the legend: https://stackoverflow.com/questions/39948151/in-r-plotly-subplot-graph-how-to-show-only-one-legend
library(plotly)
fig1 <- plot_ly(data = all_pkgs_total_downloads_TOP20,
               y = ~reorder(package, download_count), 
               x = ~download_count,
               colors = c('#d62728','#2ca02c','#9467bd','#17becf' ),
               color = ~Color_coding2,
               #split = ~Color_coding2,
               #split = ~repo,                           # color coding
               type = 'bar', orientation = 'h',
               showlegend = T,
               legendgroup=~Color_coding2
)
fig1 <- fig1 %>% layout(barmode = 'stack',
                      yaxis = list(title = ""),
                      xaxis = list(title ="downloads count"))
#fig1

fig2 <- plot_ly(data = all_pkgs_subset1_downloads_TOP10,
                y = ~reorder(package, download_count), 
                x = ~download_count,
                colors = c('#d62728','#17becf' ),
                color = ~Color_coding2,
                #split = ~Color_coding2,
                #split = ~repo,                           # color coding
                type = 'bar', orientation = 'h',
                showlegend = FALSE,
                legendgroup=~Color_coding2
              )
fig2 <- fig2 %>% layout(barmode = 'stack',
                             yaxis = list(title = ""),
                             xaxis = list(title ="downloads count"))

fig3 <- plot_ly(data = all_pkgs_subset2_downloads_TOP10,
                y = ~reorder(package, download_count), 
                x = ~download_count,
                colors = c('#d62728','#9467bd'),
                color = ~Color_coding2,
                #split = ~Color_coding2,
                #split = ~repo,                           # color coding
                type = 'bar', orientation = 'h',
                showlegend = FALSE,
                legendgroup=~Color_coding2
)

fig3 <- fig3 %>% layout(barmode = 'stack',
                            yaxis = list(title = ""),
                            xaxis = list(title ="downloads count"))

## recursive subplot as descripbed here: https://plotly-r.com/arranging-views
s1 <- subplot(fig1, nrows = 1)
s2 <- subplot(fig2, fig3, nrows = 1)



## give a title to each subplot. See thread: https://stackoverflow.com/questions/77825289/adding-titles-on-plotly-subplots-in-r
fig <- subplot(s1 , s2, nrows = 2, heights = c(0.65, 0.35),
               margin = 0.04)
fig <- fig %>% add_annotations(
  x = c(.25, .75),
  #y = 1,
  y = 0.31,
  xref = "paper",
  yref = "paper",
  text = c("(B) PCA-specific packages, top 10", "(C) Biplot-specific packages, top 10"),
  showarrow = F,
  xanchor = "center",
  yanchor = "bottom",
  font = list(size = 12)
)

fig <- fig %>% add_annotations(
  x = c(.5),
  y = 1,
  xref = "paper",
  yref = "paper",
  text = c("(A) All surveyed packages, top 20"),
  showarrow = F,
  xanchor = "center",
  yanchor = "bottom",
  font = list(size = 12)
)

## finally position and sort the legend items
## therad:https://plotly.com/r/legend/

fig <- fig %>% layout(legend = list(orientation = 'h'))

## display
fig


## [NOT WORKING] Save higher resolution - plotly only ####

## official plotly stuff DOESNT WORK
## thread https://plotly.com/r/static-image-export/
## plus a patch: https://stackoverflow.com/questions/73604954/error-when-using-python-kaleido-from-r-to-convert-plotly-graph-to-static-image
# library(reticulate)   # for saving high res
# 
# reticulate::py_install('kaleido', pip = TRUE)
# reticulate::py_require(c("kaleido"))
# # path_to_python <- C:\Users\es679\AppData\Local\r-miniconda"
# # use_python(path_to_python)
# 
# save_image(fig,
# #orca(fig,
#            paste0(Out_path,"/biplot1_highres_2.png",collapse = " ")
#            #width = 1300,
#            #height = 800
# )


## try this instead
# if ( !require(RSelenium) ) {
#   install.packages("RSelenium", repos = "https://cloud.r-project.org/")
# }
#fig %>%
#  export(file = paste0(Out_path,"/biplot1_highres_2.svg",collapse = " "),
#         selenium = RSelenium::rsDriver(browser = "chrome"))




## Save file ####
write.csv(downloads_all, paste0(Out_path,"/pkgs_downloads_mySelect.csv",collapse = " ")) 
write.csv(total_downloads, paste0(Out_path,"/pkgs_downloads_mySelect_TOT.csv",collapse = " ")) 
write.csv(all_other_downloads, paste0(Out_path,"/pkgs_downloads_everythingCRAN.csv",collapse = " "))
write.csv(all_pkgs_total_downloads, paste0(Out_path,"/pkgs_downloads_everything_TOT.csv",collapse = " ")) 