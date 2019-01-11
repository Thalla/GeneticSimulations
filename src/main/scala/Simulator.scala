import java.io.File

import Earth._
import PrintElem.PrintElem

import scala.collection.mutable.ListBuffer

/**
  * main method
  * controls biological simulation using class Cell
  * receives simulation data from Cell and transfers it to SimulationData
  * @version 4.2 jar with param, consistent param order, easier mass simulation
  */
object Simulator{
  def main(args: Array[String]): Unit = {
    //var simulator = new Simulator(args(0))

    var basePath =  "C:\\Users\\feroc\\Documents\\ThesisLocal\\SimulationTree\\"
    var simulator = new Simulator(basePath)
  }
}

class Simulator(basePath:String) {
  val steps = 100             // 1.000000 52s maxSize 100.000, 53s maxSize 50.000, 50s maxSize 500.000, 48s maxSize 1000000, 49s maxSize 2000000, 48s maxSize 1500000

    System.gc()

    //var basePath =  "C:\\Users\\feroc\\Documents\\ThesisLocal\\SimulationTree\\"
    var mrnaSeed = 2000
    var codonNumb = 64
    var geneLength = 3
    var geneNumb = 22
    var mrnaId = 0
    var initAaNumb = 20
    var similarAars = false
    var aarsSeed = 200
    var aarsLifeticksStartValue = 10
    var livingAarsSeed = 20
    var outputSeed = 1
    //var simulator = new Simulator(basePath)

    runSimulation(basePath, mrnaSeed, codonNumb, geneLength, geneNumb, mrnaId, initAaNumb, similarAars, aarsSeed, aarsLifeticksStartValue, livingAarsSeed, outputSeed)
    //Seq.range(2,20).par.foreach{changingParam => simulator.runSimulation(basePath, mrnaSeed, codonNumb, geneLength, geneNumb, mrnaId, changingParam, similarAars, aarsSeed, aarsLifeticksStartValue, livingAarsSeed, outputSeed)}

    aarsSeed = 201



  def runSimulation(basePath:String, mrnaSeed:Int, codonNumb:Int, geneLength:Int, geneNumb:Int, mrnaId:Int, initAaNumb:Int, similarAars:Boolean = true, aarsSeed:Int, aarsLifeticksStartValue:Int, livingAarsSeed:Int,  outputSeed:Int): Unit ={
   def time[R](block: => R): R = {
      val t0 = System.nanoTime()
      val result = block // call-by-name
      val t1 = System.nanoTime()
      println("Elapsed time: " + (t1 - t0) + "ns")
      result
    }

    // init Cell
    val r = new scala.util.Random(outputSeed)
    val cell = new Cell(r)
    cell.init(basePath, mrnaSeed, codonNumb, geneLength, geneNumb, mrnaId, initAaNumb, similarAars, aarsSeed, aarsLifeticksStartValue, livingAarsSeed, outputSeed)

    // init SimulationData
    val simulationData = new SimulationData(cell.path, steps)
    simulationData.init()
    val toPrint = List(PrintElem.codeTableFitness, PrintElem.aaNumb)

    // simulate
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