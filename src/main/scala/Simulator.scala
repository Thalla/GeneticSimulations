import Earth._
import PrintElem.PrintElem

import scala.collection.mutable.ListBuffer

/**
  * main method
  * controls biological simulation using class Cell
  * receives simulation data from Cell and transfers it to SimulationData
  */
object Simulator {
  val steps = 1000000         // 1.000000 52s maxSize 100.000, 53s maxSize 50.000, 50s maxSize 500.000, 48s maxSize 1000000, 49s maxSize 2000000, 48s maxSize 1500000

  def main(args: Array[String]): Unit = {

    val l:List[Double] = List.fill(20)(3.31413542346243625475251535426443765475687735623454365768)
    time(sum0)
    time(sum1())
    time(sum3())

    def sum0 ():Unit={      // 10.000000 1.829.193.554ns
      var i = 0
      while(i <= 10000000){
        l.sum
        i += 1
      }
    }


    def sum1 ():Unit={      // 10.000000 510.270.584ns
      var i = 0
      while(i <= 10000000) {
        var sum = 0.0
        val lit = l.iterator
        while (lit.hasNext) {
          sum += lit.next()
        }
        i += 1
      }
    }

    def sum3 ():Unit={      // 10.000000 10.640.276.912ns
      var i = 0
      while(i <= 10000000) {
        l.par.foldLeft(0.0){_+_}
        i += 1
      }
    }

    def time[R](block: => R): R = {
      val t0 = System.nanoTime()
      val result = block // call-by-name
      val t1 = System.nanoTime()
      println("Elapsed time: " + (t1 - t0) + "ns")
      result
    }
    val path = "C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\csv\\"
    val toPrint = List(PrintElem.generationFitness, PrintElem.aaNumb)
    val numbPrintElem = toPrint.length
    var cellData:scala.collection.mutable.Map[PrintElem, Any] = scala.collection.mutable.Map[PrintElem, Any]()
    val cell = new Cell()
    cellData += (PrintElem.generationFitness -> cell.generationFitness)
    cellData += (PrintElem.aaNumb -> cell.aaNumb)
    time(cell.init(path)) // generation -1

    time(do{
      cell.translate() // results in a new generation
      //SimulationData.updateAaHasTransl(cell.aaTranslData._2, cell.generationID)


      SimulationData.update(cell.generationFitness, cell.aaTranslData._1, cell.generationID)        ///TODO calculate fitness not while translation but while codeTable creation /((newCell.unambiguousness.foldLeft(0.0)(_+_)) /codonNumb.toDouble

      //mRNAdata = newCell.mRNAdata.reverse :: mRNAdata //10 sec per 500000 (Energiesparmodus)

    }while (cell.generationID <= steps - 2))

    SimulationData.finishOutput(List(PrintElem.generationFitness, PrintElem.aaNumb), cell.generationID)

  }
}
/*
if (cell.generationID % 10000 == 0) {
  //SimulationData.writeArrayToFile(SimulationData.generationFitness.mkString("\n"), path+"generationFitness.csv", true)

  //SimulationData.updateProtocol(cell.toHtmlString(List(PrintElem.mRNA, PrintElem.livingAARSs, PrintElem.codeTable)))
}*/