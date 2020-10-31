"This script is the main file that creates a Dash app.

Usage: app.R
"

library(tidyverse) 
library(ggplot2)
library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(dashTable)
library(plotly)
library(maps)
library(reshape2)

#Read data--------------------------------------------
W5 <- W5_data <- read_csv("W5_OTU_distribution.csv")
W1 <- read_csv("W1_Sample_Coordinates.csv")
W1$Sample <- as.factor(W1$Sample)
world_map <- map_data("world")

refKey <- read_csv("PR2_V9_IDs.csv") # Read refKey here because larger

##Keys and labels------------------------------------------

refKey <- as_tibble(refKey)

depthKey <- tibble(label=c("Surface","Deep Chlorophyll Maximum","Both"), 
                   value=c("SUR", "DCM", "Both"))


#Functions--------------------------------------------  

make_map <-function(depth='DCM') {
  
  #Get labels
  #ref_var <- "EU011928.1.1772_U"   #for testing 
  depth_var = depthKey$label[depthKey$value==depth]
  
  Data <-  W5 %>% filter(str_detect(refs, "EU011928.1.1772_U")) %>% filter(pid >= 90)
  Data <- Data[,2:335] #Get only TARA_
  Data_long <- melt(Data, variable.name="Sample", value.name="Abundance")
  Data <- left_join(Data_long, W1, by="Sample")  %>% mutate(Present = case_when(Abundance > 0 ~ "Yes", TRUE ~ "No"))
  
  
  if (depth_var=="Surface") {  #This is the label

    Data <- Data %>% filter(Depth=="SUR")
    Map <- ggplot(world_map, aes(x = long, y = lat, group = group)) +
      geom_polygon(fill="gray75", colour = "white") +
      theme_bw() +
      geom_point(Data, mapping=aes(x=long,lat, group=1, size=ifelse(Abundance==0, 0.01, Abundance), color=Present, shape=Present)) +
      scale_color_manual(values=c("gray20","red")) +
      scale_shape_manual(values=c(4, 19))
    

  } else if (depth_var=="Deep Chlorophyll Maximum")  {

    Data <- Data %>% filter(Depth=="DCM")
    Map <- ggplot(world_map, aes(x = long, y = lat, group = group)) +
      geom_polygon(fill="gray75", colour = "white") +
      theme_bw() +
      geom_point(Data, mapping=aes(x=long,lat, group=1, size=ifelse(Abundance==0, 0.01, Abundance), color=Present, shape=Present)) +
      scale_color_manual(values=c("gray20","red")) +
      scale_shape_manual(values=c(4, 19))
    
    
  } else  {
  
  
  Map <- ggplot(world_map, aes(x = long, y = lat, group = group)) +
    geom_polygon(fill="gray75", colour = "white") +
    theme_bw() +
    geom_point(Data, mapping=aes(x=long,lat, group=1, size=ifelse(Abundance==0, 0.01, Abundance), color=Present, shape=Present)) +
    scale_color_manual(values=c("gray20","red")) +
    scale_shape_manual(values=c(4, 19))
  
  }
  
  ggplotly(Map)
  
}


##Assign components of dashboard to variables---------------------


#Graphs
map <- dccGraph(
  id = "Map",
  figure = make_map()
)

#Headings and labels

heading <- htmlH1("Tara Oceans")

subtitle <- htmlH3("V9 Metabarcode Distribution")

input_label <- htmlLabel("Enter your PR2 V9 sequence ID:")

button_label <- htmlLabel("Select Depth:")


#elements

Input <- dccInput(
  id = 'Ref Input',
  type="text",
  placeholder="Enter your ID") #Need to map out options??


Button <- dccRadioItems(
  id = 'Depth Button',
  options = list(
    list("label" = "Surface", "value"="SUR"),
    list("label" = "Deep Chlorophyll Maximum", "value"="DCM"),
    list("label" = "Both", "value" = "Both")
  ),
  value='DCM'
  )


div_header <- htmlDiv(
  list(heading, 
       subtitle
  ),
  style = list(
    backgroundColor = '#46718F',
    textAlign = 'center',
    color = 'white',
    margin = 0,
    marginTop = 0,
    width= '100%'
  )
)


div_sidebar <- htmlDiv(
  list(input_label,
       Input,
       htmlBr(),
       htmlBr(),
       htmlBr(),
       button_label,
       Button
  ),
  style = list('background-color' = '#9DB7C9',
               'fontColor' = 'white',
               'padding' = 10,
               'fontSize' = 13,
               'width' = '11%')
)


div_main <- htmlDiv(
  list(map)
)

#create dash instance

app <- Dash$new()

##Dash layout----------------------------------

app$layout(
  div_header, htmlDiv(
    list(div_sidebar, div_main),
    style=list('display' = 'flex',
               'justify-content'='center',
               'width'='100%')
  )
)


##Callbacks------------------------------------

app$callback(
  output=list(id='Map',property='figure'),
  params=list(input(id='Depth Button', property='value')),
  function(depth) {
    make_map(depth)
  }
)




##Run the App---------------------------------

app$run_server(host = "0.0.0.0", ports = Sys.getenv('PORT', 8050))






