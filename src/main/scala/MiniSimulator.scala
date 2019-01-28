
object MiniSimulator{
  def main(args: Array[String]): Unit = {
    val simulator = new Simulator(args(0), args(1).toInt, args(2).toInt, args(3).toInt, args(4).toString)
    val cells = simulator.initSimulation(2000, 3, Seq(1,2,3,4,5), 1, Seq(false, true), Seq(1,2,3,4,5), Seq(1,2,3,4,5,6,7,8), Seq(1,2,3,4,5,6,7,8,9,10), 200, Seq(1,2,3,4,5), 1, true)
    simulator.runSimulation(cells)
  }
}
