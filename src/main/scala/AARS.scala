import AA._
import Earth.{getCodons}

class AARS(var aaSeq: Vector[AA], var translations:Map[(AA, Int),List[(Double, Int)]]) {

  val lifeticksStartValue = 10
  var lifeticks = 0

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

  def translationsToHtmlString():String= {
    val codons = getCodons(64)
    var content = "<table>"
    for (
      translation <- translations
    ) {
      content += "<tr>\n<td>" + translation._1._1.toString + "</td><td>" + codons(translation._1._2) + "</td><td><table>"
      for (
        t <- 0 until translation._2.length
      ) {
        content += "<tr><td>" + translation._2(t)._1 + "</td><td>" +codons(translation._2(t)._2) + "</td></tr>"
      }
      content += "</table>\n</td></tr>\n"
    }
    content += "</table>"
    content
  }
}

