import java.io.{BufferedWriter, File, FileWriter}
import scala.math.pow

/**
  * Runs a mostly predefined set of simulations
  * @version 5.0
  */
object Simulator {
  /**
    * This main Method is also used by other main objects
    * @param args path, codonNumb, livingAarsSeed, steps, translationMethod
    */
  def main(args: Array[String]): Unit = {
    if (args.length != 0) {
      val simulator = new Simulator(args(0), args(1).toInt, args(2).toInt, args(3).toInt, args(4))
      simulator.simulate()
    } else {
      val basePath = "SimulationResults\\"
      time({
        val simulator = new Simulator(basePath, 64, 20, 1000, "affinity")
        val cells: List[Map[String, Any]] = simulator.initSimulation(200, 3, Seq(20,40), Seq(false), Seq(20), Seq(6), Seq(20), 30, Seq(10), 1, true)
        simulator.runSimulation(cells, simulator.simulateCell)
      })
    }
  }

  /**
    * Timer
    * @param block
    * @tparam R
    * @return R
    */
  def time[R](block: => R): R = {
    val t0 = System.nanoTime()
    val result = block
    println("Elapsed time: " + (System.nanoTime - t0) + "ns")
    result
  }
}

/**
  * Controls biological simulation of the aaRS feedback loop. Uses objects of class Cell.
  * @param basePath path were the folders with the results of the simulation are created
  * @param codonNumb number of codons that is used for mRNA and aaRSs translations creation
  * @param livingAarsSeed seed for random number generator that decides on which living aaRSs the simulation starts with
  * @param steps number of generations that shall be simulated
  * @param translationMethod defines how the decision which translation will be taken is made in the case of an unambiguous codon
  */
class Simulator(basePath: String, codonNumb: Int, livingAarsSeed: Int, steps: Int, translationMethod: String) {
  /**
    * initializes and starts simulation
    */
  def simulate() {
    val cells: List[Map[String, Any]] = initSimulation()
    runSimulation(cells, simulateCell)
  }

  /**
    * Initializes the folder system and the following files: config, mRNA, aaRS, livingAars.
    * This is a pre-step for cell initialisation that happens in parallel afterwards.
    *
    * @param mrnaSeed seed value for the random number generator that chooses the codons used for the mRNA after all codons have been used once
    * @param geneLength defines length of aaRS sequence and is always three
    * @param geneNumbs number of genes that are translated in each cell generation
    * @param similarAarss true -> similar aaRSs sequences also have similar translations
    * @param initAaNumbs number of amino acids that are available in the simulation, maximal 23, for more amino acids Enum AA must be extended
    * @param maxAnticodonNumbs maximal number of translations an aaRSs can have
    * @param aarsLifeticksStartValues defines how long aaRSs live. An aaRS with ten lifeticks survives ten generations.
    * @param aarsSeed seed value for random number generator that defines the translations each aaRS has
    * @param livingAarsStartNumbs number of living aaRSs in the beginning of the simulation
    * @param outputSeed seed for random number generator that influences which translation is taken in the case of translation method "random" or "affinity"
    * @param addToOutput true -> the last set of living aaRSs is saved in a file. A second run of the same parameter combinations starts with the last saved set of living aaRSs.
    * @return cells: The initialization information for each cell.
    */
  def initSimulation(mrnaSeed: Int = 2000, geneLength: Int = 3, geneNumbs: Seq[Int] = Seq(2, 5, 15, 20, 30), similarAarss: Seq[Boolean] = Seq(true, false), initAaNumbs: Seq[Int] = Seq(3, 4, 8, 10, 16, 20), maxAnticodonNumbs: Seq[Int] = Seq(2, 4, 6, 8), aarsLifeticksStartValues: Seq[Int] = Seq(2, 5, 10, 15), aarsSeed: Int = 200, livingAarsStartNumbs: Seq[Int] = Seq(2, 5, 15, 20, 30), outputSeed: Int = 1, addToOutput: Boolean = true): List[Map[String, Any]] = {
    val r = new scala.util.Random(outputSeed)
    val cell = new Cell(r)
    var cells: List[Map[String, Any]] = List()

    implicit def bool2int(b: Boolean): Int = if (b) 1 else 0

    // for each parameter combination the parameter information is packed into the array "cellInfo".
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
        println("Living aaRSs start numb " + livingAarsStartNumb + " hasn't been used because it is greater than number of existing aaRSs.")
      }
    }))))))
    cells
  }

  /**
    * Runs simulations in parallel with the given simulation method
    * @param cells
    * @param simulationMethod
    */
  def runSimulation(cells: List[Map[String, Any]], simulationMethod: Map[String, Any] => Unit): Unit = {
    System.gc()
    cells.par.foreach(cellInfo => simulationMethod(cellInfo))
  }

  /**
    *
    * @param cellInfo
    */
  val simulateCell = (cellInfo: Map[String, Any]) => {
    //init cell
    val outputSeed = cellInfo("outputSeed").asInstanceOf[Int]
    val r = new scala.util.Random(outputSeed)
    val cell = new Cell(r)
    val (runFlag, outputPath) = cell.initProperties(cellInfo("translationMethod").toString, cellInfo("codonNumb").asInstanceOf[Int], cellInfo("initAaNumb").asInstanceOf[Int], cellInfo("mrnaPath").asInstanceOf[String], cellInfo("aarsPath").asInstanceOf[String], cellInfo("aarsLifeticksStartValue").asInstanceOf[Int], cellInfo("livingAarsPath").asInstanceOf[String], cellInfo("livingAarsDir").asInstanceOf[String], outputSeed, cellInfo("addToOutput").asInstanceOf[Boolean])

    if (runFlag) {
      // init SimulationData
      cell.simulationData = new SimulationData(outputPath, steps, List(PrintElem.codeTableFitness, PrintElem.livingAars))
      //cell.simulationData.updateProtocol(cell.toHtmlString(List(PrintElem.codons, PrintElem.mRNA, PrintElem.allAars, PrintElem.livingAars)))

      //translation
      do {
        cell.translationStep() // results in a new generation
        cell.simulationData.updateCodeTableFitness(cell.codeTableFitness, cell.generationID)
        //cell.simulationData.updateAaNumb(cell.aaTranslData._1, cell.generationID)
        //cell.simulationData.updateAaHasTransl(cell.aaTranslData._2, cell.generationID)
        //cell.simulationData.updateProtocol(cell.toHtmlString(List(PrintElem.livingAars, PrintElem.codeTable)))
        //mRNAdata = newCell.mRNAdata.reverse :: mRNAdata //+10 sec per 500000 (in power save mode)
      } while (cell.generationID < (steps - 1))
      //finish output
      cell.simulationData.livingAars = cell.livingAars
      cell.simulationData.finishOutput(cell.generationID)
    }
  }

  /**
    * @TODO
    */
  val simulateCellWithEvolutionaryAlg = (cellInfo: Map[String, Any]) => {

  }



}
