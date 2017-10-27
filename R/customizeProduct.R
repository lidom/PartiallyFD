customizeProduct <-
function(filename){
  if(filename == "SRL_POS_HT") {
    colSpe <- c( color = "YlGnBu", rev= FALSE)
  } else if(filename == "SRL_POS_NT") {
    colSpe <- c( color = "YlGnBu", rev= FALSE)
  } else if(filename == "SRL_NEG_HT") {
    colSpe <- c( color = "YlOrRd", rev= FALSE)
  } else{
    colSpe <- c( color = "YlOrRd", rev= FALSE)
  }
  return(list(color=colSpe))
}
