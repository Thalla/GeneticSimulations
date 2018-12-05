import java.io.{BufferedWriter, File, FileWriter}
import PrintElem._

/**
  * Receives data from Simulator
  *
  * produces file output
  */
object SimulationData {
  val path = "C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\csv\\"
  val maxSize = 1500000
  private[this] var _generationFitness: Array[Double] = new Array[Double](maxSize)
  private[this] var protocol:List[String] = List()
  private[this] var oldComparison = new Array[Array[Boolean]](20)   //TODO
  val aaChanges:Array[Array[Int]] = new Array[Array[Int]](Simulator.steps)   //genIds, added amino acids, removed amino acids
  aaChanges(0) = new Array[Int](2)
  val tableFieldOverTime:Array[Array[Int]] = new Array[Array[Int]](Simulator.steps)//Array.fill(steps)(Array.fill(8)(0))  //assumption: max numb of aaRS per field is 8 TODO
  private[this] var _aaHasTranslation: Array[Array[Boolean]] = new Array[Array[Boolean]](maxSize)
  private [this] val _aaNumb = new Array[Integer](maxSize)


  def update(generationFitnessValue:Double, aaNumbValue:Int, genID:Int):Unit = {
    if(genID%maxSize == 0 && genID != 0){
      var append = true
      if(genID == maxSize) append = false
      writeToFile(_generationFitness.mkString("\n"), path+"generationFitness.csv", append)
    }
    _generationFitness(genID%maxSize) = generationFitnessValue
  }


  /**
    *
    * @param value
    * @param genID
    */
  def updateAaNumb(value: Int, genID:Int): Unit = {
    if(genID%maxSize == 0 && genID != 0) {
      var append = true
      if(genID == maxSize) append = false
      writeToFile(_aaNumb.mkString("\n"), path+"aaNumb.csv", append)
    }
    _aaNumb(genID%maxSize) = value
  }

  /**
    *
    * @param value
    * @param genID
    */
  def updateAaHasTransl(value: Array[Boolean], genID:Int):Unit = {
    if(genID%maxSize == 0) {
      var append = true
      if(genID == 0) append = false
      writeToFile(_aaNumb.mkString("\n"), path+"aaNumb.csv", append)
    }
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
    protocol = html :: protocol
  }

  def toString(value:Array[Double], length:Int):String ={
    val strB = new StringBuilder(length)
    " "
  }

  /**
    * @param toPrint List of Tuples with filename and file content
*/
  def finishOutput(toPrint:List[PrintElem], genID:Int):Unit={

    val toFile = (data:String, filename:String) => {
      val file = new File(path + filename + ".csv")
      val bw = new BufferedWriter(new FileWriter(file, true))
      bw.newLine()
      bw.write(data)
      bw.close()
    }

    for(
      elem <- toPrint
    ){
      elem match {
        case PrintElem.generationFitness => toFile(_generationFitness.take(genID).mkString("\n"), elem.toString())
        case PrintElem.aaNumb => toFile(_aaNumb.take(genID).mkString("\n"), elem.toString())
      }
    }
  }

  def writeToFile(data:String, path:String, append:Boolean): Unit ={
    val file = new File(path)
    val bw = new BufferedWriter(new FileWriter(file, append))
    if(append == true) bw.newLine()
    bw.write(data)
    bw.close()
  }


}