
app <- Dash$new()

app$layout(htmlDiv(list(htmlH2('Hello World'),
                        dccDropdown(id = 'dropdown',
                                    options = list(
                                      list('label' = 'LA', 'value' = 'LA'),
                                      list('label' = 'NYC', 'value' = 'NYC'),
                                      list('label' = 'MTL', 'value' = 'MTL')
                                    ),
                                    value = 'LA'),
                        htmlDiv(id = 'display-value'))
)
)

app$callback(output=list(id='display-value', property='children'),
             params=list(
               input(id='dropdown', property='value')),
             function(value)
             {
               sprintf('You have selected %s', value)
             }
)


##Run the App---------------------------------

app$run_server(host = "0.0.0.0", ports = Sys.getenv('PORT', 8050))






