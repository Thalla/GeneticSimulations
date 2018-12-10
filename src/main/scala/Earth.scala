import AA._

import scala.util.Random
import java.io.{BufferedWriter, File, FileWriter}

import com.github.tototoshi.csv._

import scala.collection.mutable.{ArrayBuffer, ListBuffer}  //TODO?: replace by own methods (data-> List-> needed type => data-> needed type)


/** Fourth draft
  *
  * @nextSteps: similarity, probabilities, NN?, DM?, evo. Alg.?, analyse, parameters, Empirie, Fitness, Metrik
  * @author Hanna Schumacher
  * @version 4.1   -> 500.000 in 315s ; in 2?? s; in 272s without some stuff, 190s without all
  * TODO choose by translation fitness, not always the best one, use pathvariable
  */
object Earth extends Cell{
  //TODO: check if start parameters are valid
  //TODO: tests, exceptions, logging, class diagram, cases, directly connect code/model decisions to research
  /** initiates and starts simulation
    * @param args
    */







  /*
   * simulation; holds global parameters; cares for cell reproduction
   */
  def simulate(cell:Cell, steps:Int):Array[Array[Int]] = {



    /*











    file = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\csv\\generationFitness.csv")
    var writer = CSVWriter.open(file)
    writer.writeAll(List(generationFitness))
    writer.close()

    file = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\csv\\mRNAdata.csv")
    writer = CSVWriter.open(file)
    writer.writeAll(mRNAdata)
    writer.close()

    file = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\csv\\testdata.csv")
    bw = new BufferedWriter(new FileWriter(file))
    bw.write(aaChanges.map(_.mkString(", ")).mkString("\n"))
    bw.close()
*/


/*
    //var fitnessCSV:List[List[Double]] = List()

    /*for(
      generation:Array[Double] <- aaNumb
    ){
      fitnessCSV = fitnessCSV :+ generation.toList
    }*/
    file = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\csv\\aaNumb.csv")
    var writer = CSVWriter.open(file)
    writer.writeAll(List(aaNumbPerGen.toList))
    writer.close()







    //aaChanges to file
    var content1 =""
    var content2 =""
                            /*for(
                              i <- 0 until steps
                            ){
                              content1 += i + " "
                              content2 += i + " "
                            }
                            content1 += "\n"
                            content2 += "\n"*/
    /*for(
      i <- 0 until steps
    ){
      content1 += aaChanges(i)(0).toString() + ", "
      //TODO content2 += tableFieldOverTime(i)(0).toString() + ", "
    }
    content1 = content1.dropRight(2)
    content1 += "\n"
    for(
      i <- 0 until steps
    ){
      content1 += aaChanges(i)(1).toString() + ", "
    }
    content1 = content1.dropRight(2)
    content2 = content2.dropRight(2)
    content2 += "\n"
    for(
      i <- 0 until steps
    ){
      //TODO content2 += tableFieldOverTime(i)(1).toString() + ", "
    }
    content2 = content2.dropRight(2)
    content2 += "\n"
    for(
      i <- 0 until steps
    ){
      //TODO content2 += tableFieldOverTime(i)(2).toString() + ", "
    }
    content2 = content2.dropRight(2)
    content2 += "\n"
    for(
      i <- 0 until steps
    ){
      //TODO content2 += tableFieldOverTime(i)(3).toString() + ", "
    }
    content2 = content2.dropRight(2)
    content2 += "\n"
    for(
      i <- 0 until steps
    ){
      //TODO content2 += tableFieldOverTime(i)(4).toString() + ", "
    }
    content2 = content2.dropRight(2)
    val file1 = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\csv\\testdata.csv")
    val bw1 = new BufferedWriter(new FileWriter(file1))
    bw1.write(content1)
    bw1.close()
    val file2 = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\Modeling\\csv\\testdataFirstCell.csv")
    val bw2 = new BufferedWriter(new FileWriter(file2))
    bw2.write(content2)
    bw2.close()*/
*/
    //aaChanges
    Array()
  }




  def writeToFile(filename:String, content:String): Unit ={
    val file = new File(s"C:\\Users\\feroc\\OneDrive\\Dokumente\\HS\\Semester\\4\\Thesis\\$filename.html")
    val bw = new BufferedWriter(new FileWriter(file))
    bw.write(content)
    bw.close()
  }










  def printMRNA(mRNA:List[List[Int]], codons:IndexedSeq[(Any)]):Unit = {
    print(s"mRNA: List of genes. Gene is a List of codon IDs.")
    mRNA.foreach(gene => {
      println()
      gene.foreach(codonID => print(codons(codonID)))
    })
    println()
    println()
  }



  def print2DimTRNAmatrix (matrix:Array[Array[Any]]):Unit = {
    for (
      i <- 0 until matrix.length;
      j <- 0 until matrix(0).length
    ) {
      val elem = matrix(i)(j)
      if (elem == null) print(" ")//i + ", " + j + " ")
      else println(matrix(i)(j).toString())
    }
  }




  def listToString(list:List[Any]):String ={
    var listAsString = ""
    list.foreach(elem => listAsString += elem + " ,")
    listAsString
  }

  // prints cell list
    def toString(cells: Vector[Cell]): String ={
      var toPrint:String = ""
      for(cell <- cells) toPrint += cell.toString()
      toPrint
    }

}