import java.io.{File, IOException}

import com.github.tototoshi.csv.CSVReader
import org.junit._
import org.scalatest.junit.AssertionsForJUnit

import scala.reflect.io.Path
import scala.util.{Failure, Success, Try}

class BasicTest extends AssertionsForJUnit {




  @Before def init():Unit={
    //cell.init(path, 5, 3, 2, 16, 48, 3, 9,48, 0,false, 0)
  // Problem: ich verwende die IDs der Aminosäuren, Ala hat ID 16, aber im allAARS Array hat Ala Position 0, wenn Ala an nullter Stelle im initAA Vektor steht.

    //path = "C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\testCSVs\\mini1\\"
    /*new File(path).mkdirs
    val simulator = new Simulator(path, 9)
    simulator.initCell(path, 0, 3, 3, 2, 1, 3, false, 1, 3, 10, 2, 3, 4, true)*/

    // my minature example
    //val simulator = new Simulator(path, 3, 3, 15)
    //val cells = simulator.initSimulation(0, 3, Seq(2), 1, Seq(false), Seq(3), Seq(3), Seq(10), 1, Seq(2), 4, true)




    //29.  12:30   DI, Do, Fr
  }
  @Test def basics ():Unit = {
    var reader:CSVReader = null
    try {
      var basePath = "C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\testCSVs\\mini2\\"
      val r = new scala.util.Random(0)
      val cell = new Cell(r)

      assert(Try(Path(basePath).deleteRecursively()) match {
        case Success(lines) => true
        case Failure(f) => false
      })

      var simulator = new Simulator(basePath, 3, 3, 15)
      var cells = simulator.initSimulation(1, 3, Seq(2), 1, Seq(false), Seq(3), Seq(3), Seq(10), 1, Seq(2), 4, true)
      simulator.runSimulation(cells)

      var path: String = basePath + "mRNA_s1_c3_gl3_gn2_#1\\initAaNumb_3\\similar_false\\aaRS_s1_ac3_lt10\\livingAars_n2_s3\\output_s4\\codeTableFitness.csv"
      var file = new File(path)
      assert(file.exists()== true)
      reader = CSVReader.open(file)
      val data1 = reader.all()
      reader.close()

      assert(Try(Path(basePath).deleteRecursively()) match {
        case Success(lines) =>       true
        case Failure(f) => false
      })

      simulator = new Simulator(basePath, 3, 3, 5)
      cells = simulator.initSimulation(1, 3, Seq(2), 1, Seq(false), Seq(3), Seq(3), Seq(10), 1, Seq(2), 4, true)
      simulator.runSimulation(cells)
      simulator.runSimulation(cells)
      simulator.runSimulation(cells)
      reader = CSVReader.open(file)
      val data2 = reader.all()

      assert(data1.length == data2.length)
      assert(data1.eq(data2))
    }catch{
      case io: IOException => println("IO :(")
      case e: Throwable => println(e)
    }
    finally{
      reader.close()
    }


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
