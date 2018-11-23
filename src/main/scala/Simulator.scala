import Earth._


import scala.collection.mutable.ListBuffer

object Simulator {
  val steps = 100000

  def main(args: Array[String]): Unit = {

    var livingAARSs:ListBuffer[AARS] = Earth.init()
    val cell = new Cell(livingAARSs)

    def time[R](block: => R): R = {
      val t0 = System.nanoTime()
      val result = block // call-by-name
      val t1 = System.nanoTime()
      println("Elapsed time: " + (t1 - t0) + "ns")
      result
    }

    time(while (cell.generationID <= steps - 2) {
      val data:CellData = cell.translate()

      SimulationData.generationFitness(cell.generationID-1)= (data.numbAaWithTransl.toDouble/aaNumb.toDouble)* data.unambiguousness   ///TODO tell cell its fitness or calculate fitness not while translation but while codeTable creation /((newCell.unambiguousness.foldLeft(0.0)(_+_)) /codonNumb.toDouble
      //mRNAdata = newCell.mRNAdata.reverse :: mRNAdata //10 sec per 500000 (Energiesparmodus)

      if (cell.generationID % 10000 == 0) {
        SimulationData.addToProtocol(cell.toHtmlString(List(PrintElem.mRNA, PrintElem.livingAARSs, PrintElem.codeTable)))
      }
    })

    //SimulationData.writeToFile(List())
  }
}
