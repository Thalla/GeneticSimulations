object CustomSimulator {
  def main(args: Array[String]): Unit = {
    val simulator = new Simulator(args(0), args(1).toInt, args(2).toInt, args(3).toInt)
    val cells = simulator.initSimulation(args(4).toInt, 3, Seq(args(5).toInt), args(6).toInt, Seq(args(7).toBoolean), Seq(args(8).toInt), Seq(args(9).toInt), Seq(args(10).toInt), args(11).toInt, Seq(args(12).toInt), args(13).toInt, args(14).toBoolean)
    simulator.runSimulation(cells)
  }
}
