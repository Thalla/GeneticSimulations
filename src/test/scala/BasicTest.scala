import java.io.File

import org.junit._
import org.scalatest.junit.AssertionsForJUnit

class BasicTest extends AssertionsForJUnit {

  var path = "C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\testCSVs\\"
  val r = new scala.util.Random(0)
  val cell = new Cell(r)
  @Before def init():Unit={
    //cell.init(path, 5, 3, 2, 16, 48, 3, 9,48, 0,false, 0)
  // Problem: ich verwende die IDs der Aminosäuren, Ala hat ID 16, aber im allAARS Array hat Ala Position 0, wenn Ala an nullter Stelle im initAA Vektor steht.

    path = "C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\testCSVs\\mini1\\"
    /*new File(path).mkdirs
    val simulator = new Simulator(path, 9)
    simulator.initCell(path, 0, 3, 3, 2, 1, 3, false, 1, 3, 10, 2, 3, 4, true)*/
    val simulator = new Simulator(path, 3, 3, 15)
    val cells = simulator.initSimulation(0, 3, Seq(2), 1, Seq(false), Seq(3), Seq(3), Seq(10), 1, Seq(2), 4, true)
    simulator.runSimulation(cells)

    //29.  12:30   DI, Do, Fr
  }
  @Test def basics ():Unit = {

    /*assert(cell.aaNumb == 2)
    cell.translate()
    assert(cell.allAars(1)(0)(0).lifeticks == 10)
    assert(cell.allAars(0)(1)(0).lifeticks == 9)
    assert(cell.codeTableFitness == 0.6875)
    cell.allAars.foreach(x => x.foreach(y => y.foreach(z => println(z.toString()))))
    println(cell.mRNA.toString())
    cell.livingAars.foreach(x => println(x.toString()))*/
  }

}
