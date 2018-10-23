import AA._
import Earth.{getCodons}

class AARS(val aaSeq: Vector[AA], var translations:Map[(AA, Int),(Double, List[Int])]) {

  val lifeticksStartValue = 2
  var lifeticks = lifeticksStartValue

  def reduceLifeTicks():AARS ={
    lifeticks -= 1
    this
  }

  def resetLifeticks():Unit ={
    lifeticks = lifeticksStartValue
  }

  override def toString(): String ={
    "\n" + "aaRS: " + aaSeq.mkString(", ") + "\nlifeticks: " + lifeticks + "\nTranslations: " + translations.mkString(", ")
  }

  def translationsToHtmlString():String={
    val codons = getCodons(64)
    var content = "<table>"
    for(
      translation <- translations
    ){
      content += "<tr>\n<td>" + translation._1._1.toString + "</td><td>" + codons(translation._1._2)+"</td><td>" + codons(translation._2._2.head) + "</td></tr>\n"
    }
    content += "</table>"
    content
  }

}

