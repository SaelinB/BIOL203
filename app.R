
library(tidyverse) 
library(ggplot2)
library(dash)
library(dashCoreComponents)
library(dashHtmlComponents)
library(plotly)
library(maps)
library(dashTable)
library(reshape2)


#Read data--------------------------------------------
W5  <- read_csv("W5_OTU_Distributions_filtered2.csv")
W1 <- read_csv("W1_Sample_Coordinates.csv")
total_ab_euk <- read_csv("total_ab_euk.csv")
ID_seqs_genus <- read_csv("ID_seqs2.csv") 
W1$Sample <- as.factor(W1$Sample)
world_map <- map_data("world")


seqs_ID <- ID_seqs_genus %>% select(-Genus) %>% as.tibble() #For label-key table (as seqs - IDs)

depthKey <- tibble(label=c("Surface","Deep Chlorophyll Maximum"), 
                   value=c("SUR", "DCM"))


make_map <-function(depth='DCM', seq='blank', PID="50") {  
  

  depth_var = depthKey$label[depthKey$value==depth]
  
  IDs = seqs_ID$PR2_ID[seqs_ID$Sequence==seq]  #Get list of IDs, separated by ";"
  
  ID_var <- str_replace_all(IDs, ";", "|")  #Will become "or" statement in str_detect()
  
  PID_var <- as.numeric(PID)
  # #Make new df with totalabundances
  
  Data <-  W5 %>% filter(str_detect(refs, ID_var)) %>% filter(pid >= PID_var) #Get only TARA_ and  #Look for refs that have the ID associated with that sequence
  Data <- Data[,2:335]  #Get only TARA_
  
  Data_long <- melt(Data, variable.name="Sample", value.name="Abundance")  
  
  Data  <- left_join(Data_long, W1, by="Sample") %>%  left_join(total_ab_euk, by="Sample") %>% mutate(lat=round(lat,-1), long=round(long,-1)) %>% select(-Diversity, -size_range, -Sample)
  
  #Sample is multiple of # OTUs = summing OTU abundance, averaging total reads for each (rather than summing total multiple times) 
  Data_tot_euk  <- Data %>% group_by(lat, long, Depth) %>% summarise("Sum_ab" = sum(Abundance), "Total_Reads" = sum(unique(Totab_all))) %>%
    mutate(per_tot_reads = Sum_ab/Total_Reads) %>% mutate(per_tot_read_ad=round(per_tot_reads*1000000, 2))
  
  
  
  if (ID_var =="blank") {    #need to do this to keep dots small, otherwise are very large when size=Abundance or ifelse
    
    Map <- ggplot(world_map, aes(x = long, y = lat, group = group)) +
      geom_polygon(fill="gray75", colour = "white") +   #world map
      theme_bw() +
      geom_point(Data_tot_euk, mapping=aes(x=long, lat, group=1), color="grey90") +
      labs(x="Longitude", y="Latitude")
    
    ggplotly(Map, width=1500, height=900) 
    
  } else if (depth_var=="Surface") {  #This is the label
    
    
    Data <- Data_tot_euk %>% filter(Sum_ab > 0) %>% filter(Depth=="SUR") 
    
    if (nrow(Data) == 0) {
      
      Map <- ggplot(world_map, aes(x = long, y = lat, group = group)) +
        geom_polygon(fill="gray75", colour = "white") +   #world map
        theme_bw() +
        labs(x="Longitude", y="Latitude")
      
      ggplotly(Map, width=1500, height=900,  tooltip = c("text")) %>%
        layout(title=list(text="No Hits", y=0.99), font=list(size=12))
      
    } else {
      
      Map <- ggplot(world_map, aes(x = long, y = lat, group = group)) +
        geom_polygon(fill="gray75", colour = "white") +  #world map
        theme_bw() +
        geom_point(Data, mapping=aes(x=long, lat,  group=1, size=per_tot_reads, 
                                     text = paste("Long:", long, "<br>", "Lat:", lat, "<br>",  "Percent Total Reads", per_tot_reads, "<br>", "Number of Reads:", Sum_ab)), alpha=0.6, color="red") +
        labs(x="Longitude", y="Latitude") 
      
      ggplotly(Map, width=1500, height=900,  tooltip = c("text")) %>% 
        add_annotations(x = Data$long, y = Data$lat,showarrow =FALSE, text = Data$per_tot_read_ad, clicktoshow="onoff", yshift=14) %>%
        layout(title=list(text="OTU Abundance (percent total eukaryotic barcodes x10<sup>-6</sup>)", y=0.99), font=list(size=12))
      
    }
    
    
  }  else  {
    
    Data <- Data_tot_euk %>% filter(Sum_ab > 0) %>% filter(Depth=="DCM") 
    
    
    if (nrow(Data) == 0) {
      
      Map <- ggplot(world_map, aes(x = long, y = lat, group = group)) +
        geom_polygon(fill="gray75", colour = "white") +   #world map
        theme_bw() +
        labs(x="Longitude", y="Latitude")
      
      ggplotly(Map, width=1500, height=900,  tooltip = c("text")) %>%
        layout(title=list(text="No Hits", y=0.99), font=list(size=12))
      
    } else {
      
      
      Map <- ggplot(world_map, aes(x = long, y = lat, group = group)) +
        geom_polygon(fill="gray75", colour = "white") +  #world map
        theme_bw() +
        geom_point(Data, mapping=aes(x=long, lat,  group=1, size=per_tot_reads,
                                     text = paste("Long:", long, "<br>", "Lat:", lat, "<br>",  "Percent Total Reads", per_tot_reads, "<br>", "Number of Reads:", Sum_ab)), alpha=0.6, color="red") +
        labs(x="Longitude", y="Latitude") 
      
      ggplotly(Map, width=1500, height=900,  tooltip = c("text")) %>% 
        add_annotations(x = Data$long, y = Data$lat,showarrow =FALSE, text = Data$per_tot_read_ad, clicktoshow="onoff", yshift=14) %>%
        layout(title=list(text="OTU Abundance (percent total eukaryotic barcodes x10<sup>-6</sup>)", y=0.99), font=list(size=12))
      
    }
  
   }
    
}





##Assign components of dashboard to variables---------------------

#Graphs
map <- dccGraph(
  id = "Map",
  figure = make_map()
)

#Headings and labels

heading <- htmlH2("Tara Oceans V9 Metabarcode Distribution", style=list("font-family"= "arial"))

input_label <- htmlLabel("Enter your V9 sequence:", style=list("font-family"= "arial"))

PIDinput_label <- htmlLabel("Minimum Percent Identity to OTU:", style=list("font-family"="arial"))

button_label <- htmlLabel("Select Depth:", style=list("font-family"="arial"))

#elements

Input <- dccInput(
  id = 'Seq Input',
  type="text",
  placeholder="") 

pidInput <- dccInput(
  id = 'PID Input',
  type="text",
  placeholder="") 


Button <- dccRadioItems(
  id = 'Depth Button',
  options = list(
    list("label" = "Surface", "value"="SUR"),
    list("label" = "Deep Chlorophyll Maximum", "value"="DCM")
  ), style=list("font-family"= "arial"),
  value='DCM'
)


div_header <- htmlDiv(
  list(heading),
  style = list(
    backgroundColor ='white',
    textAlign = 'center',
    color = 'black')
)


div_sidebar <- htmlDiv(
  list(input_label,
       htmlBr(),
       Input,  #Seq input
       htmlBr(),
       htmlBr(),
       PIDinput_label,
       htmlBr(),
       pidInput, #PID input
       htmlBr(),
       htmlBr(),
       button_label,
       Button,
       htmlBr()
  ),
  style = list(backgroundColor ='white',
               color = 'black',
               'padding' = 10,
               'fontSize' = 15)
)


div_main <- htmlDiv(
  list(map))



#create dash instance
app <- Dash$new()

##Dash layout----------------------------------

app$layout(
  htmlDiv(list(div_header, div_sidebar,div_main))
)




##Callbacks------------------------------------


app$callback(
  output=list(id='Map',property='figure'),
  params=list(input(id='Depth Button', property='value'),
              input(id='Seq Input', property='value'),
              input(id='PID Input', property='value')),
  function(depth, seq, PID) {
    make_map(depth, seq, PID)
  }
)



##Run the App---------------------------------
app$run_server(host = "0.0.0.0", ports = Sys.getenv('PORT', 8050))


