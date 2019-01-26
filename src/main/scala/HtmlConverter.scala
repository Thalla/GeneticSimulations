object HtmlConverter {
  val break = "</br>"
  val newLine = "\n"

  val header = (number:Int, text:Any)=>{
    s"<h$number> ${text.toString()} </h$number>$newLine"
  }

  val p = (text:Any)=>{
    s"<p>${text.toString()}</p>$newLine"
  }

  val divColour = (text:Any, numb:Int, numbAll:Int)=>{
    val colour = (255/4)*numb
    //val colourNegative = 255-colour
    s"<div style='font-weight: bold; display: inline-block; background-color: rgb(${colour%255}, ${colour*20%255}, ${colour*300%255}); color: rgb(${255-colour%255}, ${255-colour*20%255}, ${255-colour*300%255})'>${text.toString}</div>"
  }

  val listElement = (text:Any)=>{
    s"<li>$newLine${text.toString}</li>$newLine"
  }

  val list = (text:Any)=>{
    s"<ul>$newLine${text.toString()}</ul>$newLine"
  }

  val table = (text:Any)=>{
    s"<table>$newLine${text.toString}</table>$newLine"
  }

  val tableHeader = (text:Any)=>{
    s"<th>$newLine${text.toString}</th>$newLine"
  }

  val tableRow = (text:Any)=>{
    s"<tr>$newLine${text.toString}</tr>$newLine"
  }

  def tableField(text:Any):String={
    s"<td>$newLine${text.toString}</td>$newLine"
  }
}
