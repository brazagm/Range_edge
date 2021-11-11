# function made by Eran Raviv 2015
# for facilitate addition of textbox in plot
# located by the corners

Corner_text <- function(text, location="topright"){
  legend(location,legend=text, bty ="n", pch=NA) 
}