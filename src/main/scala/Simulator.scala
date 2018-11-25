import Earth._


import scala.collection.mutable.ListBuffer

object Simulator {
  val steps = 500000

  def main(args: Array[String]): Unit = {



    def time[R](block: => R): R = {
      val t0 = System.nanoTime()
      val result = block // call-by-name
      val t1 = System.nanoTime()
      println("Elapsed time: " + (t1 - t0) + "ns")
      result
    }
    val path = "C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\csv\\"
    val cell = new Cell()
    time(cell.init(path))
    val celly = new Cell()
    time(celly.init(path))
    /*time(while (cell.generationID <= steps - 2) {
      cell.translate()
      //SimulationData.aaTranslData_=(cell.aaTranslData)
      //SimulationData.generationFitness(cell.generationID-1)= (data.numbAaWithTransl.toDouble/aaNumb.toDouble)* data.unambiguousness   ///TODO tell cell its fitness or calculate fitness not while translation but while codeTable creation /((newCell.unambiguousness.foldLeft(0.0)(_+_)) /codonNumb.toDouble
      //mRNAdata = newCell.mRNAdata.reverse :: mRNAdata //10 sec per 500000 (Energiesparmodus)

      if (cell.generationID % 10000 == 0) {
        //SimulationData.addToProtocol(cell.toHtmlString(List(PrintElem.mRNA, PrintElem.livingAARSs, PrintElem.codeTable)))
      }
    })*/

    var i = 0
    time(while(i <= 5000000){
      cell.translate()
      i += 1
    })
    i = 0
    time(while(i <= 5000000){
      celly.yTranslate()
      i += 1
    })

    //SimulationData.writeToFile(List())
  }
}
