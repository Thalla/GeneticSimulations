import AA._
import HtmlConverter._

class AARS(var aaSeq: Vector[AA], var translations:Map[(AA, Int),List[(Double, Int)]], val lifeticksStartValue:Int) {

  var lifeticks = 0

  def reduceLifeTicks():AARS ={
    lifeticks -= 1
    this
  }

  def resetLifeticks():Unit ={
    lifeticks = lifeticksStartValue
  }

  override def toString(): String ={
    var aaSeqIt = aaSeq.iterator
    var aaSeqIds:List[Int] = List()
      while(aaSeqIt.hasNext){
        val aa:AA = aaSeqIt.next()
        aaSeqIds = aa.id +: aaSeqIds
      }
    s"${aaSeqIds.reverse.mkString(",")},$lifeticks,${translations.mkString(",")}"
  }

  /*override def toString(): String ={
    "\n" + "aaRS: " + aaSeq.mkString(", ") + "\nlifeticks: " + lifeticks + "\nTranslations: " + translations.mkString(", ")
  }*/

  def translationsToHtmlString(codons:Array[Any]):String= {
    val codonNumb = codons.length
    var tableContent = ""
    for (
      translation <- translations
    ) {
      tableContent += tableRow(tableField(translation._1._1.toString) + tableField(divColour(codons(translation._1._2), translation._1._2, codonNumb)) + tableField(translation._2(0)._1) + tableField(codons(translation._2(0)._2)))
      for (
        t <- 1 until translation._2.length
      ) {
        tableContent += tableRow(tableField("") + tableField("") + tableField(translation._2(t)._1) + tableField(codons(translation._2(t)._2)))
      }
    }
     table(tableContent)
  }
}

