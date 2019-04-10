object MiniSimulator{
  /**
    * Most parameters are predefined with low values
    * @param args basePath, codonNumb, livingAarsSeed, steps, translationMethod
    */
  def main(args: Array[String]): Unit = {
    val simulator = new Simulator(basePath = args(0), codonNumb = args(1).toInt, livingAarsSeed = args(2).toInt, steps = args(3).toInt, translationMethod = args(4).toString)
    val cells = simulator.initSimulation(geneNumbs = Seq(1,2,3,4,5), similarAarss = Seq(false, true), initAaNumbs = Seq(1,2,3,4,5), maxAnticodonNumbs = Seq(1,2,3,4,5,6,7,8), aarsLifeticksStartValues = Seq(1,2,3,4,5,6,7,8,9,10), livingAarsStartNumbs = Seq(1,2,3,4,5))
    simulator.runSimulation(cells, simulator.simulateCell)
  }
}
