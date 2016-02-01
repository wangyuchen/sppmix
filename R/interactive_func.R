#' @export
Get_User_Input_sppmix<- function(prompt_string="",modeYN=1)
{
  options(warn=-1)
  #modeYN=1, ask for yes or no
  #modeYN=0, ask for a value
  if(modeYN)
  {
    ret=0;#no, 1 is yes
    while(1)
    {
      ANSWER <- readline(paste(prompt_string,
                               "(Y)es or (N)o? or (Q)uit "))
      if (substr(ANSWER, 1, 1) == "n"
          || substr(ANSWER, 1, 1) == "N")
      {
        ret=0
        break
      }
      if (substr(ANSWER, 1, 1) == "y"
          || substr(ANSWER, 1, 1) == "Y")
      {
        ret=1
        break
      }
      if (substr(ANSWER, 1, 1) == "q"
          || substr(ANSWER, 1, 1) == "Q")
      {
        stop("Execution ended by the user")
      }
    }
  }
  else
  {
    ret=0;#default return value
    while(1)
    {
      #message("Enter ",prompt_string,":")
      val <- readline(paste("Enter",prompt_string,": "))
      #scan(what=double())
      #      check=is.na(as.numeric(val))
      if(is.na(as.numeric(val)))
        message("Enter a number, not letters")
      else
      {
        ret=as.numeric(val)
        break
      }
    }
  }
  options(warn=0)
  return(ret)
}
