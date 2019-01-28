object CustomSimulator {
  def main(args: Array[String]): Unit = {
    val simulator = new Simulator(args(0), args(1).toInt, args(2).toInt, args(3).toInt, args(4).toString)
    val cells = simulator.initSimulation(args(5).toInt, 3, Seq(args(6).toInt), args(7).toInt, Seq(args(8).toBoolean), Seq(args(9).toInt), Seq(args(10).toInt), Seq(args(11).toInt), args(12).toInt, Seq(args(13).toInt), args(14).toInt, args(15).toBoolean)
    simulator.runSimulation(cells)
  }
}
