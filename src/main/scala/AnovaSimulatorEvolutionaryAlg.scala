object AnovaSimulatorEvolutionaryAlg {
  def main(args: Array[String]): Unit = {

    if (args.length != 0) {
      val simulator = new Simulator(basePath = args(0), codonNumb = args(1).toInt, livingAarsSeed = args(2).toInt, steps = args(3).toInt, translationMethod = args(4).toString)
      val cells = simulator.initSimulation(geneNumbs = Seq(5, 10, 20, 30), similarAarss = Seq(false, true), initAaNumbs = Seq(args(5).toInt), maxAnticodonNumbs = Seq(6), aarsLifeticksStartValues = Seq(2, 5, 10, 15), livingAarsStartNumbs = Seq(5, 10, 20, 30), outputSeed = args(6).toInt)
      simulator.runSimulation(cells, simulator.simulateCellWithEvolutionaryAlg)

    } else {
      val basePath = "C:\\Users\\feroc\\Documents\\ThesisLocal\\SimulationTree3\\translMethod_affinity\\"
      val params: Array[String] = Array(basePath, "64", "20", "11", "affinity", "23", "3")
      AnovaSimulatorEvolutionaryAlg.main(params)
    }
  }
}

