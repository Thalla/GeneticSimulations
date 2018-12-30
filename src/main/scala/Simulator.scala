import java.io.File

import Earth._
import PrintElem.PrintElem

import scala.collection.mutable.ListBuffer

/**
  * main method
  * controls biological simulation using class Cell
  * receives simulation data from Cell and transfers it to SimulationData
  */
object Simulator {
  val steps = 1000000             // 1.000000 52s maxSize 100.000, 53s maxSize 50.000, 50s maxSize 500.000, 48s maxSize 1000000, 49s maxSize 2000000, 48s maxSize 1500000

  def main(args: Array[String]): Unit = {

    System.gc()

    Seq.range(0,10).par.foreach{ i => runSimulation(i)}



  }

  def runSimulation(i:Int): Unit ={
    val path = s"C:\\Users\\feroc\\Documents\\ThesisLocal\\Simulations\\$i\\"
    val f = new File(path).mkdirs()
    val simulationData = new SimulationData(path)

    def time[R](block: => R): R = {
      val t0 = System.nanoTime()
      val result = block // call-by-name
      val t1 = System.nanoTime()
      println("Elapsed time: " + (t1 - t0) + "ns")
      result
    }

    val toPrint = List(PrintElem.codeTableFitness, PrintElem.aaNumb)
    val numbPrintElem = toPrint.length
    var cellData:scala.collection.mutable.Map[PrintElem, Any] = scala.collection.mutable.Map[PrintElem, Any]()
    val cell = new Cell()
    cellData += (PrintElem.codeTableFitness -> cell.codeTableFitness)
    cellData += (PrintElem.aaNumb -> cell.aaNumb)
    time(cell.init(path, newAARS = true, newLivingAARS = true, similarAARS = false, initAA = Vector(AA.Phe, AA.Leu, AA.Ser, AA.Tyr, AA.Cys, AA.Trp), newMRNA = true)) // generation -1
    simulationData.init()

    time(do{
      cell.translate() // results in a new generation
      //SimulationData.updateAaHasTransl(cell.aaTranslData._2, cell.generationID)


      simulationData.update(cell.codeTableFitness, cell.aaTranslData._1, cell.generationID)        ///TODO calculate fitness not while translation but while codeTable creation /((newCell.unambiguousness.foldLeft(0.0)(_+_)) /codonNumb.toDouble

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