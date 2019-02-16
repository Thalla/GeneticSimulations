import java.io.{BufferedWriter, File, FileWriter}

import PrintElem._

import scala.collection.mutable.ListBuffer

/**
  * Receives data from Simulator
  *
  * produces file output
  */
class SimulationData(path:String, steps:Int, toPrint:List[PrintElem]) {
  val maxSize = 1500000
  private[this] var _codeTableFitness: Array[Double] = new Array[Double](maxSize)
  private[this] var protocol:List[String] = List()
  private[this] var oldComparison = new Array[Array[Boolean]](20)   //TODO
  val aaChanges:Array[Array[Int]] = new Array[Array[Int]](steps)   //genIds, added amino acids, removed amino acids
  aaChanges(0) = new Array[Int](2)
  val tableFieldOverTime:Array[Array[Int]] = new Array[Array[Int]](steps)//Array.fill(steps)(Array.fill(8)(0))  //assumption: max numb of aaRS per field is 8 TODO
  private[this] var _aaHasTranslation: Array[Array[Boolean]] = new Array[Array[Boolean]](maxSize)
  private [this] val _aaNumb = new Array[Integer](maxSize)
  var livingAars:ListBuffer[AARS] = ListBuffer()

  init()


  def init(): Unit ={

      if(!new File(path + "codeTableFitness.csv").exists()){
        for(
          elem <- toPrint
        ) {
          elem match {
            case PrintElem.codeTableFitness => writeToFile("", path+"codeTableFitness.csv", false, false)
            case PrintElem.aaNumb => writeToFile("", path+"aaNumb.csv", false, false)
            case _ =>
          }
        }
      }
      else {
        for(
          elem <- toPrint
        ) {
          elem match {
            case PrintElem.codeTableFitness => writeToFile("", path + "codeTableFitness.csv", true, true)
            case PrintElem.aaNumb => writeToFile("", path + "aaNumb.csv", true, true)
            case _ =>
          }
        }
      }
  }

  /**
    *
    * @param value new fitness value to be added
    * @param genID update position
    */
  def updateCodeTableFitness(value:Double, genID:Int):Unit = {
    if(genID%maxSize == 0 && genID != 0){
      var append = true
      if(genID == maxSize) append = false
      writeToFile(_codeTableFitness.mkString("\n"), path+"codeTableFitness.csv", append, append)
    }
    _codeTableFitness(genID%maxSize) = value
  }


  /**
    *
    * @param value new aaNumb value to be added
    * @param genID update position
    */
  def updateAaNumb(value: Int, genID:Int): Unit = {
    if(genID%maxSize == 0 && genID != 0) {
      var append = true
      if(genID == maxSize) append = false
      writeToFile(_aaNumb.mkString("\n"), path+"aaNumb.csv", append, append)
    }
    _aaNumb(genID%maxSize) = value
  }

  /**
    *
    * @param value new aaNumb value to be added
    * @param genID update position
    */
  def updateAaHasTransl(value: Array[Boolean], genID:Int):Unit = {
    if(genID%maxSize == 0 && genID != 0) {
      var append = true
      if(genID == maxSize) append = false
      writeToFile(_aaHasTranslation.mkString("\n"), path+"aaHasTranslation.csv", append, append)
    }
    _aaHasTranslation(genID%maxSize) = value
  }

  /**
    *
    * @param fieldPos
    * @param codeTable
    * @param genID
    */
  def updateTableFieldObserver(fieldPos:Tuple2[Int, Int], codeTable:Array[Array[List[AARS]]], genID:Int):Unit= {
    // if field not empty
    if (codeTable(fieldPos._1)(fieldPos._2) != null) {
      val fieldElemNumb = codeTable(fieldPos._1)(fieldPos._2).length  //numb of AARS in table field
      val field = codeTable(fieldPos._1)(fieldPos._2)
      tableFieldOverTime(genID) = new Array[Int](fieldElemNumb)
      for(
        aarsPos <- 0 until fieldElemNumb  //foreach aaRS
      ){
        tableFieldOverTime(genID)(fieldElemNumb)= field(aarsPos).aaSeq(0).id + field(aarsPos).aaSeq(1).id * 20 + field(aarsPos).aaSeq(2).id * 40    //value is calculated aaRS id
      }
    }
  }

  /**
    *
    * @param html
    */
  def updateProtocol(html:String): Unit ={
    //protocol = html :: protocol
    writeToFile(html, path + "protocol.html", true, true)
  }

  def toString(value:Array[Double], length:Int):String ={
    val strB = new StringBuilder(length)
    " "
  }

  def finishOutput(genID:Int):Unit={
    for(
      elem <- toPrint
    ){
      elem match {
        case PrintElem.codeTableFitness => writeToFile(_codeTableFitness.take(genID+1).mkString("\n"), path + elem.toString() + ".csv", true, false)
        case PrintElem.aaNumb => writeToFile(_aaNumb.take(genID+1).mkString("\n"), path + elem.toString() + ".csv", true, false)
        case PrintElem.livingAars => {
          var i = 0   // exists in any case because there are already files in this path therefore a set of living aaRS as well.
          while ((new File(path + s"livingAars${i}.csv").exists())) {
            i += 1
          }
          writeToFile(livingAars.mkString("\n"), path + s"${elem.toString}$i.csv", false, false)
        }
        case PrintElem.protocol => writeToFile(protocol.mkString("\n"), path + elem.toString + ".html", true, false)
      }
    }
  }

  def writeToFile(data:String, path:String, append:Boolean, newLine:Boolean): Unit ={
    val file = new File(path)
    val bw = new BufferedWriter(new FileWriter(file, append))
    bw.write(data)
    if(newLine) bw.newLine()
    bw.close()
  }


}
