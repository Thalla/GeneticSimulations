/**
  * main method
  * controls biological simulation using class Cell
  * receives simulation data from Cell and transfers it to SimulationData
  * @version 4.2 jar with param, consistent param order, easier mass simulation
  */
object Simulator{

  def main(args: Array[String]): Unit = {
    new Simulator(args(0), args(1).toInt)

    //val basePath =  "C:\\Users\\feroc\\Documents\\ThesisLocal\\SimulationTree1\\"
    //new Simulator(basePath, 20)
  }
}

class Simulator(basePath:String, livingAarsSeed:Int) {
  val steps = 100001             // 1.000000 52s maxSize 100.000, 53s maxSize 50.000, 50s maxSize 500.000, 48s maxSize 1000000, 49s maxSize 2000000, 48s maxSize 1500000

    System.gc()

      //var basePath =  "C:\\Users\\feroc\\Documents\\ThesisLocal\\SimulationTree\\"
    var mrnaSeed = 2000
      var codonNumb = 64
    var geneLength = 3
      //var geneNumb = 22
    var mrnaId = 1
     /*var initAaNumb = 20
     var similarAars = false*/
    var aarsSeed = 200
      /*var maxAnticodonNumb = 6
      var aarsLifeticksStartValue = 10
      var livingAarsStartNumb = 22*/
        //var livingAarsSeed = 20
    var outputSeed = 1
    var addToOutput = true

    var cells:Seq[Cell] = Seq()


    //initCell(basePath, mrnaSeed, codonNumb, geneLength, geneNumb, mrnaId, initAaNumb, similarAars, aarsSeed, maxAnticodonNumb, aarsLifeticksStartValue, livingAarsStartNumb, livingAarsSeed, outputSeed, addToOutput)
    //Seq.range(2,20).par.foreach{changingParam => simulator.runSimulation(basePath, mrnaSeed, codonNumb, geneLength, geneNumb, mrnaId, changingParam, similarAars, aarsSeed, maxAnticodonNumb, aarsLifeticksStartValue, livingAarsStartNumb, livingAarsSeed, outputSeed, addToOutput)}

    //Vector(48, 64).foreach(codonNumb => Seq(2, 5, 15, 20, 30).foreach(geneNumb => Seq(3, 4, 8, 10, 16, 20).foreach(initAaNumb => Seq(true, false).foreach(similarAars => Seq(2, 4, 6, 8).foreach(maxAnticodonNumb => Seq(2, 5, 10, 15).foreach(aarsLifeticksStartValue => Seq(2, 5, 15, 20, 30).foreach(livingAarsStartNumb => initCell(basePath, mrnaSeed, codonNumb, geneLength, geneNumb, mrnaId, initAaNumb, similarAars, aarsSeed, maxAnticodonNumb, aarsLifeticksStartValue, livingAarsStartNumb, livingAarsSeed, outputSeed, addToOutput))))))))
    //Vector(48, 64).foreach(codonNumb => Seq(2, 5, 15, 20, 30).foreach(geneNumb => Seq(3, 4, 8, 10, 16, 20).foreach(initAaNumb => Seq(true, false).foreach(similarAars => Seq(2, 4, 6, 8).foreach(maxAnticodonNumb => Seq(2, 5, 10, 15).foreach(aarsLifeticksStartValue => Seq(2, 5, 15, 20, 30).foreach(livingAarsStartNumb => runSimulation(cells))))))))

  codonNumb = 48
  //
  Seq(2, 5, 15, 20, 30).foreach(geneNumb => Seq(3, 4, 8, 10, 16, 20).foreach(initAaNumb => Seq(true, false).par.foreach(similarAars => Seq(2, 4, 6, 8).foreach(maxAnticodonNumb => Seq(2, 5, 10, 15).foreach(aarsLifeticksStartValue => Seq(2, 5, 15, 20, 30).foreach(livingAarsStartNumb => {
    initCell(basePath, mrnaSeed, codonNumb, geneLength, geneNumb, mrnaId, initAaNumb, similarAars, aarsSeed, maxAnticodonNumb, aarsLifeticksStartValue, livingAarsStartNumb, livingAarsSeed, outputSeed, addToOutput)
    runSimulation(cells)
    cells = Seq()
  }))))))


  codonNumb = 64
  Seq(2, 5, 15, 20, 30).foreach(geneNumb => Seq(3, 4, 8, 10, 16, 20).foreach(initAaNumb => Seq(true, false).foreach(similarAars => Seq(2, 4, 6, 8).foreach(maxAnticodonNumb => Seq(2, 5, 10, 15).foreach(aarsLifeticksStartValue => Seq(2, 5, 15, 20, 30).foreach(livingAarsStartNumb => {
    initCell(basePath, mrnaSeed, codonNumb, geneLength, geneNumb, mrnaId, initAaNumb, similarAars, aarsSeed, maxAnticodonNumb, aarsLifeticksStartValue, livingAarsStartNumb, livingAarsSeed, outputSeed, addToOutput)
    runSimulation(cells)
    cells = Seq()
  }))))))

  codonNumb = 32
  Seq(2, 5, 15, 20, 30).foreach(geneNumb => Seq(3, 4, 8, 10, 16, 20).foreach(initAaNumb => Seq(true, false).foreach(similarAars => Seq(2, 4, 6, 8).foreach(maxAnticodonNumb => Seq(2, 5, 10, 15).foreach(aarsLifeticksStartValue => Seq(2, 5, 15, 20, 30).foreach(livingAarsStartNumb => {
    initCell(basePath, mrnaSeed, codonNumb, geneLength, geneNumb, mrnaId, initAaNumb, similarAars, aarsSeed, maxAnticodonNumb, aarsLifeticksStartValue, livingAarsStartNumb, livingAarsSeed, outputSeed, addToOutput)
    runSimulation(cells)
    cells = Seq()
  }))))))

  def runSimulation(cells:Seq[Cell]): Unit ={
    // simulate
    cells.par.foreach(cell => {do {
      cell.translate() // results in a new generation
      cell.simulationData.updateCodeTableFitness(cell.codeTableFitness, cell.generationID) //((newCell.unambiguousness.foldLeft(0.0)(_+_)) /codonNumb.toDouble
      cell.simulationData.updateAaNumb(cell.aaTranslData._1, cell.generationID)
      cell.simulationData.updateAaHasTransl(cell.aaTranslData._2, cell.generationID)
      //mRNAdata = newCell.mRNAdata.reverse :: mRNAdata //10 sec per 500000 (Energiesparmodus)
    } while (cell.generationID <= steps - 2)

    cell.simulationData.livingAars = cell.livingAars
    cell.simulationData.finishOutput(List(PrintElem.codeTableFitness, PrintElem.aaNumb, PrintElem.livingAars), cell.generationID)
      this.cells = cells.filterNot(_.equals(cell))
  })
  }


  def initCell(basePath:String, mrnaSeed:Int, codonNumb:Int, geneLength:Int, geneNumb:Int, mrnaId:Int, initAaNumb:Int, similarAars:Boolean, aarsSeed:Int, maxAnticodonNumb:Int, aarsLifeticksStartValue:Int, livingAarsStartNumb:Int, livingAarsSeed:Int,  outputSeed:Int, addToOutput:Boolean): Unit = {
    /*def time[R](block: => R): R = {
      val t0 = System.nanoTime()
      val result = block // call-by-name
      val t1 = System.nanoTime()
      println("Elapsed time: " + (t1 - t0) + "ns")
      result
    }*/

    // init Cell
    val r = new scala.util.Random(outputSeed)
    val cell = new Cell(r)
    val runFlag: Boolean = cell.init(basePath, mrnaSeed, codonNumb, geneLength, geneNumb, mrnaId, initAaNumb, similarAars, aarsSeed, maxAnticodonNumb:Int, aarsLifeticksStartValue, livingAarsStartNumb, livingAarsSeed, outputSeed, addToOutput)

    if (runFlag) {

    // init SimulationData
    cell.simulationData = new SimulationData(cell.path, steps)
    cell.simulationData.init()
    cell.toPrint = List(PrintElem.codeTableFitness, PrintElem.aaNumb)

      cells = cell +: cells



    }
  }
}




/*
if (cell.generationID % 10000 == 0) {
  //SimulationData.writeArrayToFile(SimulationData.codeTableFitness.mkString("\n"), path+"codeTableFitness.csv", true)

  //SimulationData.updateProtocol(cell.toHtmlString(List(PrintElem.mRNA, PrintElem.livingAARSs, PrintElem.codeTable)))
}*/