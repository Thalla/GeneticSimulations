import java.io.{BufferedWriter, File, FileWriter}

/**
  * main method
  * controls biological simulation using class Cell
  * receives simulation data from Cell and transfers it to SimulationData
  *
  * @version 4.3 first write to file system afterwards only read
  */
object Simulator {

  def main(args: Array[String]): Unit = {
    //new Simulator(args(0), args(1).toInt)

    val basePath = "C:\\Users\\feroc\\Documents\\ThesisLocal\\SimulationTree2\\selection_false\\translMethod_probability\\"
    new Simulator(basePath, 64, 20, 11)
  }
}

class Simulator(basePath: String, codonNumb: Int, livingAarsSeed: Int, steps: Int) {

  var cells: List[Map[String, Any]] = initSimulation()
  simulate(cells)

  /**
    * Defines the file structure and initializes the following files: config, mRNA, aaRS, livingAars.
    *
    * @param mrnaSeed
    * @param geneLength
    * @param mrnaId
    * @param aarsSeed
    * @param outputSeed
    * @param addToOutput
    */
  def initSimulation(mrnaSeed: Int = 2000, geneLength: Int = 3, mrnaId: Int = 1, aarsSeed: Int = 200, outputSeed: Int = 1, addToOutput: Boolean = true): List[Map[String, Any]] = {
    val r = new scala.util.Random(outputSeed)
    val cell = new Cell(r)
    var cells: List[Map[String, Any]] = List()

    implicit def bool2int(b: Boolean): Int = if (b) 1 else 0

    Seq(2, 5, 15, 20, 30).foreach(geneNumb => Seq(3, 4, 8, 10, 16, 20).foreach(initAaNumb => Seq(true, false).foreach(similarAars => Seq(2, 4, 6, 8).foreach(maxAnticodonNumb => Seq(2, 5, 10, 15).foreach(aarsLifeticksStartValue => Seq(2, 5, 15, 20, 30).foreach(livingAarsStartNumb => {
      var path = basePath
      var newConfig = false

      val codons = cell.getCodons(codonNumb)
      val (mrnaName, newCombiM) = cell.writeNewMrna(path, mrnaSeed, codons, geneLength, geneNumb, mrnaId)
      path += s"$mrnaName\\"
      val mrnaPath = path + s"$mrnaName.csv"
      if (newCombiM == true) newConfig = true


      path += s"initAaNumb_$initAaNumb\\"

      path += s"similar_$similarAars\\"
      val (aarsName, newCombiA) = cell.writeNewAars(path, similarAars, aarsSeed, codons, geneLength, geneNumb, maxAnticodonNumb, aarsLifeticksStartValue)
      path += s"$aarsName\\"
      val aarsPath = path + s"$aarsName.csv"
      if (newCombiA == true) newConfig = true

      //init livingAARS
      val (livingAarsName, newCombiL) = cell.writeNewLivingAars(path, livingAarsSeed, livingAarsStartNumb, initAaNumb)
      path += s"$livingAarsName\\"
      val livingAarsPath = path + s"$livingAarsName.csv"
      if (newCombiL == true) newConfig = true


      val cellInfo = Map("codonNumb" -> codonNumb, "initAaNumb" -> initAaNumb, "mrnaPath" -> mrnaPath, "aarsPath" -> aarsPath, "aarsLifeticksStartValue" -> aarsLifeticksStartValue, "livingAarsPath" -> livingAarsPath, "outputSeed" -> outputSeed, "addToOutput" -> addToOutput)
      cells = cellInfo +: cells

      if (newConfig) {
        val configData: List[String] = List(mrnaSeed.toString, codonNumb.toString, geneLength.toString, geneNumb.toString, mrnaId.toString, initAaNumb.toString, (similarAars: Int).toString, maxAnticodonNumb.toString, aarsLifeticksStartValue.toString, livingAarsStartNumb.toString, livingAarsSeed.toString, outputSeed.toString, (addToOutput: Int).toString, mrnaPath, aarsPath, livingAarsPath)
        val file = new File(basePath + "config.csv")
        val bw = new BufferedWriter(new FileWriter(file, true))
        bw.write(configData.mkString(","))
        bw.newLine()
        bw.close()
      }

    }))))))
    cells
  }

  def simulate(cells: List[Map[String, Any]]): Unit = {
    System.gc()

    cells.par.foreach(cellInfo => {
      val outputSeed = cellInfo("outputSeed").asInstanceOf[Int]
      val r = new scala.util.Random(outputSeed)
      val cell = new Cell(r)
      cell.path = basePath

      val runFlag = cell.initProperties(cellInfo("codonNumb").asInstanceOf[Int], cellInfo("initAaNumb").asInstanceOf[Int], cellInfo("mrnaPath").asInstanceOf[String], cellInfo("aarsPath").asInstanceOf[String], cellInfo("aarsLifeticksStartValue").asInstanceOf[Int], cellInfo("livingAarsPath").asInstanceOf[String], outputSeed, cellInfo("addToOutput").asInstanceOf[Boolean])

      if (runFlag) {
        // init SimulationData
        cell.simulationData = new SimulationData(cell.path, steps)
        cell.simulationData.init()
        cell.toPrint = List(PrintElem.codeTableFitness, PrintElem.aaNumb)

        {
          cell.translate() // results in a new generation
          cell.simulationData.updateCodeTableFitness(cell.codeTableFitness, cell.generationID) //((newCell.unambiguousness.foldLeft(0.0)(_+_)) /codonNumb.toDouble
          cell.simulationData.updateAaNumb(cell.aaTranslData._1, cell.generationID)
          cell.simulationData.updateAaHasTransl(cell.aaTranslData._2, cell.generationID)
          //mRNAdata = newCell.mRNAdata.reverse :: mRNAdata //10 sec per 500000 (Energiesparmodus)
        }
        while (cell.generationID <= steps - 2)
          cell.simulationData.livingAars = cell.livingAars
          cell.simulationData.finishOutput(List(PrintElem.codeTableFitness, PrintElem.aaNumb, PrintElem.livingAars), cell.generationID)
      }
    })
  }
}
