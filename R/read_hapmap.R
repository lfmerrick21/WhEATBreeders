read_hapmap <- function(infile)
{
  output=fread(infile)
  output[output=="NA"]<-"NA"
  output$`assembly#`=as.character(output$`assembly#`)
  output$center=as.character(output$center)
  output$protLSID=as.character(output$protLSID)
  output$assayLSID=as.character(output$assayLSID)
  output$panelLSID=as.character(output$panelLSID)
  output$QCcode=as.character(output$QCcode)

  return(output)
}
