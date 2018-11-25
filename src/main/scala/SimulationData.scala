import java.io.{BufferedWriter, File, FileWriter}
import PrintElem._

object SimulationData {
  val path = "C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\Scala\\"

  var mRNAdata:List[List[Int]] = List()
  val generationFitness= new Array[Double](Simulator.steps)
  var protocol:List[String] = List()
  var oldComparison = new Array[Boolean](Earth.aaNumb)
  //oldComparison = cell.aaHasTransl
  var newComparison = new Array[Boolean](Earth.aaNumb)
  val aaChanges:Array[Array[Int]] = new Array[Array[Int]](Simulator.steps)   //genIds, added amino acids, removed amino acids
  aaChanges(0) = new Array[Int](2)
  val tableFieldOverTime:Array[Array[Int]] = new Array[Array[Int]](Simulator.steps)//Array.fill(steps)(Array.fill(8)(0))  //assumption: max numb of aaRS per field is 8 TODO
  private[this] var _aaTranslData: (Int, Array[Boolean]) = (0, Array())

  def aaTranslData_=(value: (Int, Array[Boolean])): Unit = {
    _aaTranslData = value
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
  def addToProtocol(html:String): Unit ={
    protocol = html :: protocol
  }

  def toString(value:Array[Double], length:Int):String ={
    val strB = new StringBuilder(length)
    " "
  }

  /**
    * @param toPrint List of Tuples with filename and file content
    */
  def writeToFile(toPrint:List[Tuple2[String, String]]):Unit={
    for(
      elem <- toPrint
    ){
      val file = new File(path + elem._1)
      val bw = new BufferedWriter(new FileWriter(file))
      bw.write(elem._2)
      bw.close()
    }
  }
}
