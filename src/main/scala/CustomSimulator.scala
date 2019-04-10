object CustomSimulator {
  /**
    * Simulates one parameter set.
    *
    * @param args all simulation parameters: basePath, codonNumb, livingAarsSeed, steps, translationMethod, mrnaSeed, geneNumbs, similarAarss, initAaNumbs, maxAnticodonNumb, aarsLifeticksStartValue, livingAarsStartNumbs, outputSeed, addToOutput
    */
  def main(args: Array[String]): Unit = {
    val simulator = new Simulator(basePath = args(0), codonNumb = args(1).toInt, livingAarsSeed = args(2).toInt, steps = args(3).toInt, translationMethod = args(4).toString)
    val cells = simulator.initSimulation(mrnaSeed = args(5).toInt, geneNumbs = Seq(args(6).toInt), similarAarss = Seq(args(7).toBoolean), initAaNumbs = Seq(args(8).toInt), maxAnticodonNumbs = Seq(args(9).toInt), aarsLifeticksStartValues = Seq(args(10).toInt), livingAarsStartNumbs = Seq(args(11).toInt), outputSeed = args(12).toInt, addToOutput = args(13).toBoolean)
    simulator.runSimulation(cells, simulator.simulateCell)
  }
}
