import Earth._
import PrintElem.PrintElem

import scala.collection.mutable.ListBuffer

/**
  * main method
  * controls biological simulation using class Cell
  * receives simulation data from Cell and transfers it to SimulationData
  */
object Simulator {
  val steps = 1000000         // 1.000000 52s maxSize 100.000, 53s maxSize 50.000, 50s maxSize 500.000, 48s maxSize 1000000, 49s maxSize 2000000, 48s maxSize 1500000

  def main(args: Array[String]): Unit = {



    def time[R](block: => R): R = {
      val t0 = System.nanoTime()
      val result = block // call-by-name
      val t1 = System.nanoTime()
      println("Elapsed time: " + (t1 - t0) + "ns")
      result
    }
    val path = "C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\csv\\"
    val toPrint = List(PrintElem.generationFitness, PrintElem.aaNumb)
    val numbPrintElem = toPrint.length
    var cellData:scala.collection.mutable.Map[PrintElem, Any] = scala.collection.mutable.Map[PrintElem, Any]()
    val cell = new Cell()
    cellData += (PrintElem.generationFitness -> cell.generationFitness)
    cellData += (PrintElem.aaNumb -> cell.aaNumb)
    time(cell.init(path)) // generation -1

    time(do{
      cell.translate() // results in a new generation
      //SimulationData.updateAaHasTransl(cell.aaTranslData._2, cell.generationID)


      SimulationData.update(cell.generationFitness, cell.aaTranslData._1, cell.generationID)        ///TODO calculate fitness not while translation but while codeTable creation /((newCell.unambiguousness.foldLeft(0.0)(_+_)) /codonNumb.toDouble

      //mRNAdata = newCell.mRNAdata.reverse :: mRNAdata //10 sec per 500000 (Energiesparmodus)

    }while (cell.generationID <= steps - 2))

    SimulationData.finishOutput(List(PrintElem.generationFitness, PrintElem.aaNumb), cell.generationID)

  }
}
/*
if (cell.generationID % 10000 == 0) {
  //SimulationData.writeArrayToFile(SimulationData.generationFitness.mkString("\n"), path+"generationFitness.csv", true)

  //SimulationData.updateProtocol(cell.toHtmlString(List(PrintElem.mRNA, PrintElem.livingAARSs, PrintElem.codeTable)))
}*/