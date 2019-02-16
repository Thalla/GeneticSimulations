import java.io.{BufferedWriter, File, FileWriter}
import scala.math.pow

/**
  * controls biological simulation using class Cell
  * receives simulation data from Cell and transfers it to SimulationData
  *
  * @version 4.5 failure corrections, some documentation added
  * @TODO always finish output folder even if program is stopped
  * @TODO codonNumb = 3 -> sth. else can't be higher than ... ???
  */
object Simulator {
  /**
    * Runs a mostly predefined set of simulations
    * @param args path, codonNumb, livingAarsSeed, translationMethod
    */
  def main(args: Array[String]): Unit = {
    if(args.length != 0){
      val simulator = new Simulator(args(0), args(1).toInt, args(2).toInt, args(3).toInt, args(4))
      simulator.simulate()
    }else{
      val basePath = "C:\\Users\\feroc\\Documents\\ThesisLocal\\SimulationTree3\\translMethod_affinity\\"
      val params :Array[String] = Array(basePath, "64", "20", "11", "affinity", "23", "3")
      AnovaSimulator.main(params)
      //val simulator = Simulator(basePath, 64, 20, 11, "affinity")
      //simulator.simulate()
    }


  }
}

/**
  *
  * @param basePath
  * @param codonNumb
  * @param livingAarsSeed
  * @param steps
  * @param translationMethod
  */
class Simulator(basePath: String, codonNumb: Int, livingAarsSeed: Int, steps: Int, translationMethod:String) {
  /**
    * initializes and starts simulation
    */
  def simulate() {
    val cells: List[Map[String, Any]] = initSimulation()
    runSimulation(cells)
  }

  /**
    * Defines the file structure and initializes the following files: config, mRNA, aaRS, livingAars.
    * @param mrnaSeed
    * @param geneLength
    * @param geneNumbs
    * @param similarAarss
    * @param initAaNumbs
    * @param maxAnticodonNumbs
    * @param aarsLifeticksStartValues
    * @param aarsSeed
    * @param livingAarsStartNumbs
    * @param outputSeed
    * @param addToOutput
    * @return
    */
  def initSimulation(mrnaSeed: Int = 2000, geneLength: Int = 3, geneNumbs: Seq[Int] = Seq(2, 5, 15, 20, 30), similarAarss: Seq[Boolean] = Seq(true, false), initAaNumbs: Seq[Int] = Seq(3, 4, 8, 10, 16, 20), maxAnticodonNumbs: Seq[Int] = Seq(2, 4, 6, 8), aarsLifeticksStartValues: Seq[Int] = Seq(2, 5, 10, 15), aarsSeed: Int = 200, livingAarsStartNumbs: Seq[Int] = Seq(2, 5, 15, 20, 30), outputSeed: Int = 1, addToOutput: Boolean = true): List[Map[String, Any]] = {
    val r = new scala.util.Random(outputSeed)
    val cell = new Cell(r)
    var cells: List[Map[String, Any]] = List()

    implicit def bool2int(b: Boolean): Int = if (b) 1 else 0

    geneNumbs.foreach(geneNumb => initAaNumbs.foreach(initAaNumb => similarAarss.foreach(similarAars => maxAnticodonNumbs.foreach(maxAnticodonNumb => aarsLifeticksStartValues.foreach(aarsLifeticksStartValue => livingAarsStartNumbs.foreach(livingAarsStartNumb => {
      //param validity check
      if (!(livingAarsStartNumb > pow(geneLength, initAaNumb))) {
        var path = basePath
        var newConfig = false

        val codons = cell.getCodons(codonNumb)
        val (mrnaName, newCombiM) = cell.writeNewMrna(path, mrnaSeed, codons, geneLength, geneNumb)
        path += s"$mrnaName\\"
        val mrnaPath = path + s"$mrnaName.csv"
        if (newCombiM) newConfig = true

        path += s"initAaNumb_$initAaNumb\\"

        path += s"similar_$similarAars\\"
        val (aarsName, newCombiA) = cell.writeNewAars(path, similarAars, aarsSeed, codonNumb, initAaNumb, maxAnticodonNumb, aarsLifeticksStartValue)
        path += s"$aarsName\\"
        val aarsPath = path + s"$aarsName.csv"
        if (newCombiA) newConfig = true

        //init livingAARS
        val (livingAarsName, newCombiL) = cell.writeNewLivingAars(path, livingAarsSeed, livingAarsStartNumb, initAaNumb)
        path += s"$livingAarsName\\"
        val livingAarsDir = path
        val livingAarsPath = path + s"$livingAarsName.csv"
        if (newCombiL) newConfig = true

        val cellInfo = Map("translationMethod" -> translationMethod, "codonNumb" -> codonNumb, "initAaNumb" -> initAaNumb, "mrnaPath" -> mrnaPath, "aarsPath" -> aarsPath, "aarsLifeticksStartValue" -> aarsLifeticksStartValue, "livingAarsPath" -> livingAarsPath, "livingAarsDir" -> livingAarsDir, "outputSeed" -> outputSeed, "addToOutput" -> addToOutput)
        cells = cellInfo +: cells

        val outputPath = livingAarsDir + s"\\output_s$outputSeed\\"
        if (!new File(outputPath).exists()) newConfig = true

        if (newConfig) {
          val configData: List[String] = List(translationMethod, mrnaSeed.toString, codonNumb.toString, geneLength.toString, geneNumb.toString, initAaNumb.toString, (similarAars: Int).toString, maxAnticodonNumb.toString, aarsLifeticksStartValue.toString, livingAarsStartNumb.toString, livingAarsSeed.toString, outputSeed.toString, (addToOutput: Int).toString, mrnaPath, aarsPath, livingAarsPath, livingAarsDir)
          val file = new File(basePath + "config.csv")
          val bw = new BufferedWriter(new FileWriter(file, true))
          bw.write(configData.mkString(","))
          bw.newLine()
          bw.close()
        }
      } else {
        println("Living aaRS start numb" + livingAarsStartNumb + " hasn't been used because it is greater than number of existing aaRS.")
      }
    }))))))
    cells
  }

  /**
    *
    * @param cells
    */
  def runSimulation(cells: List[Map[String, Any]]): Unit = {
    System.gc()

    cells.par.foreach(cellInfo => simulateCell(cellInfo))
  }

  /**
    *
    * @param cellInfo
    */
  def simulateCell(cellInfo: Map[String, Any]): Unit = {
    //init cell
    val outputSeed = cellInfo("outputSeed").asInstanceOf[Int]
    val r = new scala.util.Random(outputSeed)
    val cell = new Cell(r)
    cell.path = basePath
    val (runFlag, outputPath) = cell.initProperties(cellInfo("translationMethod").toString, cellInfo("codonNumb").asInstanceOf[Int], cellInfo("initAaNumb").asInstanceOf[Int], cellInfo("mrnaPath").asInstanceOf[String], cellInfo("aarsPath").asInstanceOf[String], cellInfo("aarsLifeticksStartValue").asInstanceOf[Int], cellInfo("livingAarsPath").asInstanceOf[String], cellInfo("livingAarsDir").asInstanceOf[String], outputSeed, cellInfo("addToOutput").asInstanceOf[Boolean])

    if (runFlag) {
      // init SimulationData
      cell.simulationData = new SimulationData(outputPath, steps, List(PrintElem.codeTableFitness, PrintElem.livingAars))
      //cell.simulationData.updateProtocol(cell.toHtmlString(List(PrintElem.codons, PrintElem.mRNA, PrintElem.allAars, PrintElem.livingAars)))

      //translation
      do {
        cell.translationStep() // results in a new generation
        cell.simulationData.updateCodeTableFitness(cell.codeTableFitness, cell.generationID) //((newCell.unambiguousness.foldLeft(0.0)(_+_)) /codonNumb.toDouble
        //cell.simulationData.updateAaNumb(cell.aaTranslData._1, cell.generationID)
        //cell.simulationData.updateAaHasTransl(cell.aaTranslData._2, cell.generationID)
        //cell.simulationData.updateProtocol(cell.toHtmlString(List(PrintElem.livingAars, PrintElem.codeTable)))
        //mRNAdata = newCell.mRNAdata.reverse :: mRNAdata //10 sec per 500000 (Energiesparmodus)
      } while (cell.generationID < (steps - 1))
      //finish output
      cell.simulationData.livingAars = cell.livingAars
      cell.simulationData.finishOutput(cell.generationID)
    }
  }
}
