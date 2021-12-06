# The Natural Mortality Tool (NMT)

The Natrual Mortality tool offers an accessible way to estimate natural mortality fro a variety of life history information. The tool also offers a means to include uncertainty and combine multiple estimators into a non-parametric distribution. 
<br></br>

The tool can be launched directly from this site: https://connect.fisheries.noaa.gov/natural-mortality-tool/

You can also locally use this code by downloading the code and installing the below libraries.

## Installing libraries and locally running the NMT tool
```R
packages<-c("shiny","fishmethods","ggplot2","truncnorm","data.table","RColorBrewer","viridis","reshape2")

Running the tool can be accomplished in any of the following ways:
1) shiny::runApp(ENTER HERE USER PATH TO FOLDER CONTAINING THE NMT files)
2) Open the server.r or ui.r files in RStudio and push the "Run App" button (top rigt corner of the source panel). 
	I recommend using the "Run External" option within the "Run App" button (see small arrow in button to change options)
3) runGitHub("Natural-Mortality-Tool", "shcaba",destdir=mydir) where mydir is the path you chose to obtain results.

```

