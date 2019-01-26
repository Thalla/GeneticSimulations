import Simulator._


class MiniSimulator(basePath: String, codonNumb: Int, livingAarsSeed: Int, steps: Int) {
  def main(args: Array[String]): Unit = {
    val r = new scala.util.Random(0)
    val cell = new Cell(r)
    val simulator = new Simulator(basePath, 3, 3, 15)
    val cells = simulator.initSimulation(2000, 3, Seq(1,2,3,4,5), 1, Seq(false, true), Seq(1,2,3,4,5), Seq(1,2,3,4,5,6,7,8), Seq(1,2,3,4,5,6,7,8,9,10), 200, Seq(1,2,3,4,5), 1, true)
    simulator.runSimulation(cells)
  }
}
