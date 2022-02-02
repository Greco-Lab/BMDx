library(shiny)
library(shinyjs)

runApp(
	appDir="/home/docker/bmdx", 
	port=8080, 
	launch.browser=FALSE, 
	host="0.0.0.0"
)
