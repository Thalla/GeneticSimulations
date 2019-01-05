import java.io.File

import Earth._
import PrintElem.PrintElem

import scala.collection.mutable.ListBuffer

/**
  * main method
  * controls biological simulation using class Cell
  * receives simulation data from Cell and transfers it to SimulationData
  * @version 4.1 issue: livingAARS aren't distinct in file. -> first collection, make it distinct, then to file
  */
object Simulator {
  val steps = 100             // 1.000000 52s maxSize 100.000, 53s maxSize 50.000, 50s maxSize 500.000, 48s maxSize 1000000, 49s maxSize 2000000, 48s maxSize 1500000

  def main(args: Array[String]): Unit = {

    System.gc()

    Seq.range(0,1).par.foreach{ outputSeed => runSimulation(outputSeed)}

  }

  def runSimulation(outputSeed:Int): Unit ={
   def time[R](block: => R): R = {
      val t0 = System.nanoTime()
      val result = block // call-by-name
      val t1 = System.nanoTime()
      println("Elapsed time: " + (t1 - t0) + "ns")
      result
    }

    val r = new scala.util.Random(outputSeed)
    val cell = new Cell(r)
    val basePath = "C:\\Users\\feroc\\Documents\\ThesisLocal\\SimulationTree\\"
    cell.init(basePath, 22, 0, 22, 3, 10, 3, 23, 48, 0, true, outputSeed)

    val simulationData = new SimulationData(cell.path)
    simulationData.init()
    val toPrint = List(PrintElem.codeTableFitness, PrintElem.aaNumb)

    time(do{
      cell.translate() // results in a new generation
      simulationData.updateCodeTableFitness(cell.codeTableFitness, cell.generationID)        //((newCell.unambiguousness.foldLeft(0.0)(_+_)) /codonNumb.toDouble
      simulationData.updateAaNumb(cell.aaTranslData._1, cell.generationID)
      simulationData.updateAaHasTransl(cell.aaTranslData._2, cell.generationID)
      //mRNAdata = newCell.mRNAdata.reverse :: mRNAdata //10 sec per 500000 (Energiesparmodus)
    }while (cell.generationID <= steps - 2))

    simulationData.finishOutput(List(PrintElem.codeTableFitness, PrintElem.aaNumb), cell.generationID)
  }
}




/*
if (cell.generationID % 10000 == 0) {
  //SimulationData.writeArrayToFile(SimulationData.codeTableFitness.mkString("\n"), path+"codeTableFitness.csv", true)

  //SimulationData.updateProtocol(cell.toHtmlString(List(PrintElem.mRNA, PrintElem.livingAARSs, PrintElem.codeTable)))
}*/